function [N_opt, alpha_opt, AST_opt, AIL_at_opt] = optim_non_adaptive(params, ch_type)
% OPTIM_NON_ADAPTIVE  Non-adaptive AST maximisation (Sec. V-A).
%
%   [N_opt, alpha_opt, AST_opt, AIL_at_opt] = OPTIM_NON_ADAPTIVE(params, ch_type)
%
%   Solves the non-adaptive Average Secrecy Throughput (AST) maximisation:
%
%       max_{N, α, γ_T}   T̄ = (m(1−ε)/N) · Pr{γ̃_b ≥ γ_T}
%       subject to:        δ̄(N, α, γ_T) ≤ φ       (AIL constraint)
%                          γ_T ≥ 0
%                          N ∈ [N_min, N_max]
%                          α ∈ (0, 1]
%
%   The design variable γ_T is the on-off transmission threshold: Alice
%   transmits only when her main-channel SNR exceeds γ_T.
%
%   Solver: MATLAB's fmincon (interior-point) via the optim problem framework
%
%   Inputs:
%     params   - cell array: {M, epsilon, k, SNR, mu_B, mu_E, Nmax, Nmin,
%                             phi, Kb_factor, Ke_factor}
%     ch_type  - string; 'Rayl' or 'Rice'
%
%   Outputs:
%     N_opt       - optimal blocklength N*
%     alpha_opt   - optimal power allocation α*
%     AST_opt     - optimal AST value T̄*  [bpcu]
%     AIL_at_opt  - AIL value at the optimal point (verify ≤ φ)
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Sec. V-A.
%
% =========================================================================

[M, epsilon, k, SNR, mu_B, mu_E, Nmax, Nmin, phi, Kb_factor, Ke_factor] = params{:};

% ── Warm-start: find a feasible initial γ_T at N_max, α = 1 ─────────────
x   = sym('x');
Vb  = channel_dispersion(inf);
Rh  = M/Nmax + sqrt(Vb/Nmax) * qfuncinv(epsilon);
fx  = (1 + x) / (exp(Rh * log(2))) - 1;

if strcmp(ch_type, 'Rayl')
    cdf_sym = 1 - exp(-fx / SNR / mu_E);
else
    cdf_sym = 1 - generalized_marcum_q(sqrt(2*k*Ke_factor), ...
                     sqrt(2*(1 + Ke_factor) * fx / SNR / mu_E), k);
end

gamma_init = double(vpasolve(cdf_sym == 1 - phi, x));
if isempty(gamma_init), gamma_init = 0; end
clearvars x fx cdf_sym;

% ── Optimisation problem ──────────────────────────────────────────────────
N      = optimvar('N',      'Type', 'continuous', 'LowerBound', Nmin, 'UpperBound', Nmax);
alpha  = optimvar('alpha',  'LowerBound', 0, 'UpperBound', 1);
gammaT = optimvar('gammaT', 'LowerBound', 0);

initialPoint.N      = Nmax;
initialPoint.alpha  = 1;
initialPoint.gammaT = gamma_init;

problem = optimproblem('ObjectiveSense', 'Maximize');

% Objective: AST = m(1−ε)·Pr{γ̃_b ≥ γ_T} / N
    function obj = objective(N, gammaT)
        if strcmp(ch_type, 'Rayl')
            obj = M * (1 - epsilon) * (1 - gammainc(gammaT/SNR/mu_B, k)) / N;
        else
            a = sqrt(2 * k * Kb_factor);
            b = sqrt(2 * (1 + Kb_factor) * gammaT / SNR / mu_B);
            obj = M * (1 - epsilon) * generalized_marcum_q(a, b, k) / N;
        end
    end

% Constraint: AIL ≤ φ  (plus feasibility guards)
    function cst = constraint(N, alpha, gammaT)
        Vb_inf  = channel_dispersion(inf);
        Rinf    = M/N + sqrt(Vb_inf/N) * qfuncinv(epsilon);
        x0      = (1 + alpha*gammaT) / (exp(Rinf*log(2))) - 1;
        if strcmp(ch_type, 'Rayl')
            tau = 1 + (x0*(1-alpha)) / (alpha*(k-1));
            ail = 1 - exp(-x0/alpha/SNR/mu_E) / (tau^(k-1));
            cst = [1 - ail <= 0.9*phi;   % AIL ≥ 1−0.9φ  (lower safety guard)
                   0 <= 1 - ail;          % 1−AIL ≥ 0
                   0 <= x0;              % x₀ ≥ 0
                   0 <= ail];            % AIL ≥ 0
        else
            cst = generalized_marcum_q(sqrt(2*k*Ke_factor), ...
                      sqrt(2*(1+Ke_factor)*x0/SNR/mu_E), k) <= phi;
        end
    end

problem.Objective           = fcn2optimexpr(@objective,  N, gammaT, ch_type);
problem.Constraints.cst     = constraint(N, alpha, gammaT);

options = optimoptions('fmincon', ...
    'Algorithm',            'interior-point', ...
    'EnableFeasibilityMode', true, ...
    'SubproblemAlgorithm',  'cg', ...
    'MaxIterations',        1e5, ...
    'MaxFunctionEvaluations', 1e4, ...
    'Display',              'off');

[sol, fval, ~] = solve(problem, initialPoint, 'Options', options);

if isfinite(fval) && isreal(fval)
    N_opt      = floor(sol.N);
    alpha_opt  = sol.alpha;
    AST_opt    = fval;

    % Verify AIL constraint at solution
    Vb_inf  = channel_dispersion(inf);
    R_opt   = M/N_opt + sqrt(Vb_inf/N_opt) * qfuncinv(epsilon);
    x0_opt  = (1 + alpha_opt*sol.gammaT) / (exp(R_opt*log(2))) - 1;
    tau_opt = 1 + (x0_opt*(1-alpha_opt)) / (alpha_opt*(k-1));
    AIL_at_opt = exp(-x0_opt/alpha_opt/SNR/mu_E) / (tau_opt^(k-1));

    if AIL_at_opt <= phi
        fprintf('SUCCESS: AIL = %.2e ≤ φ = %.2e,  AST* = %.4f bpcu\n', ...
            AIL_at_opt, phi, AST_opt);
    else
        fprintf('FAILED:  AIL = %.2e > φ = %.2e,  AST* = %.4f bpcu\n', ...
            AIL_at_opt, phi, AST_opt);
    end
else
    N_opt      = Nmax;
    alpha_opt  = 1;
    AST_opt    = 0;
    AIL_at_opt = NaN;
    warning('optim_non_adaptive: solver did not converge.');
end

end
