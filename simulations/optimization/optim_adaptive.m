function [N_opt, alpha_opt, AST_opt, AST_init] = optim_adaptive(params, gammahB, ch_type)
% OPTIM_ADAPTIVE  Adaptive AST maximisation using instantaneous CSI (Sec. V-B).
%
%   [N_opt, alpha_opt, AST_opt, AST_init] = OPTIM_ADAPTIVE(params, gammahB, ch_type)
%
%   For the adaptive scheme, Alice knows the instantaneous main-channel SNR
%   γ̃_b = γ̄_b‖h_b‖² and optimises (N, α) per channel realisation:
%
%       max_{N, α}   (m(1−ε)/N) · 𝟙[δ̄(N, α, γ_b) ≤ φ]
%
%   Since the objective is discontinuous (indicator function), a genetic
%   algorithm (GA) is used as the solver.
%
%   Inputs:
%     params   - cell array: {M, epsilon, SNR, mu_B, mu_E, Nmax, Nmin, k, phi}
%     gammahB  - scalar; instantaneous normalised channel gain ‖h_b‖² for Bob
%     ch_type  - string; 'Rayl' or 'Rice'
%
%   Outputs:
%     N_opt    - optimal blocklength N*  (integer)
%     alpha_opt- optimal power allocation α*
%     AST_opt  - optimal objective value T̄*  [bpcu]
%     AST_init - objective at the initial feasible point
%
%   Solver: MATLAB Genetic Algorithm (ga) via Global Optimization Toolbox
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Sec. V-B.
%
% =========================================================================

[M, epsilon, SNR, mu_B, mu_E, Nmax, Nmin, k, phi] = params{:};

% ── Optimisation variables ────────────────────────────────────────────────
N     = optimvar('N',    'Type', 'integer', 'LowerBound', Nmin, 'UpperBound', Nmax);
alpha = optimvar('alpha', 'LowerBound', 0, 'UpperBound', 1);

% Warm-start at (N_max, α=1)
initialPoint.N     = Nmax;
initialPoint.alpha = 1;

% Evaluate initial point for reference
if strcmp(ch_type, 'Rayl')
    gammaB_init = 1 * SNR * mu_B * gammahB;
    p_init = {Nmax, M, epsilon, 1, mu_E, gammaB_init, SNR, k, 1, 1};
    [ail_init, ~] = ail_approx(p_init, 'Rayl', 'ANI');
else
    ail_init = 1;
end
AST_init = (M * (1-epsilon) * (ail_init <= phi)) / Nmax;

% ── Objective function ────────────────────────────────────────────────────
    function obj = objective_fn(N, alpha)
        gammaB = alpha * SNR * mu_B * gammahB;
        Vb     = channel_dispersion(gammaB);
        R      = M/N + sqrt(Vb/N) * qfuncinv(epsilon);
        x0     = max((1 + gammaB) / (2^R) - 1, 0);

        if strcmp(ch_type, 'Rayl')
            tau    = 1 + (x0*(1-alpha)) / (alpha*(k-1));
            ail    = exp(-x0/alpha/SNR/mu_E) / (tau^(k-1));
        else
            ail = 1;   % placeholder for Rician (extend as needed)
        end

        % Discontinuous objective: transmit only if AIL ≤ φ
        obj = (M * (1-epsilon) / N) * (ail^2 <= phi^2);
    end

% ── Problem definition ────────────────────────────────────────────────────
problem = optimproblem('ObjectiveSense', 'Maximize');
problem.Objective = fcn2optimexpr(@objective_fn, N, alpha);

options = optimoptions('ga', 'ConstraintTolerance', 1e-6, 'Display', 'off');
[solution, fval] = solve(problem, initialPoint, 'Solver', 'ga', 'Options', options);

% ── Extract results ───────────────────────────────────────────────────────
if ~isempty(solution.N)
    N_opt     = solution.N;
    alpha_opt = solution.alpha;
    AST_opt   = fval;
else
    N_opt     = Nmax;
    alpha_opt = 1;
    AST_opt   = 0;
    warning('optim_adaptive: GA solver returned empty solution.');
end

end
