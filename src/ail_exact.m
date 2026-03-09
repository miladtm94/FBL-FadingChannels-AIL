function ail = ail_exact(params, ch_type, beamforming_type)
% AIL_EXACT  Exact Average Information Leakage via numerical integration.
%
%   ail = AIL_EXACT(params, ch_type, beamforming_type)
%
%   Evaluates the exact AIL δ̄ by numerically integrating eq. (10):
%
%       δ̄ = ∫_{x₀}^{∞}  Q(√(N/V_e(x)) · [log₂((1+γ_b)/(1+x)) − R₀])
%                         · f_{γ_e}(x) dx
%
%   where Q(·) is the Gaussian Q-function, V_e(x) is the channel dispersion
%   of Eve at SNR x, and f_{γ_e} is the PDF of Eve's received SNR.
%
%   Inputs:
%     params            - cell array:
%                           {N, M, epsilon, alpha, mu_E, gammaB, SNR, k,
%                            Ke_factor, modelParams}
%     ch_type           - string; 'Rayl' or 'Rice'
%     beamforming_type  - string; 'ANI' or 'MRT'
%
%   Output:
%     ail               - scalar; exact AIL δ̄ (Monte-Carlo validated)
%
%   Note:
%     For Rician fading, the exact γ_e distribution is intractable in closed
%     form (eq. 33). Instead, a fitted parametric PDF (Gamma or Rician) is
%     used as a high-accuracy proxy — validated in Fig. 2 of the paper.
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Eq. (10).
%
% =========================================================================

[N, M, epsilon, alpha, mu_E, gammaB, SNR, k, ~, modelParams] = params{:};

Rs = M / N;
Vb = channel_dispersion(gammaB);
R  = Rs + sqrt(Vb ./ N) .* qfuncinv(epsilon);   % effective codeword rate
x0 = eval_point(gammaB, 'normal', {N, epsilon, M});

% ── Integrand building blocks ─────────────────────────────────────────────
Ve       = @(x) log2(exp(1)).^2 .* (1 - (1 + x).^(-2));
func1    = @(x) sqrt(N ./ Ve(x)) .* (log2((1 + gammaB) ./ (1 + x)) - R);
myQfunc  = @(x) 0.5 * erfc(x / sqrt(2));   % Q(x) = 0.5·erfc(x/√2)

% Helper: robust integral that retries with shrinking upper limit if NaN
function val = safe_integral(f, lo, hi)
    val = integral(f, lo, hi, 'RelTol', 1e-4, 'AbsTol', 1e-6);
    if isnan(val)
        exp_limit = 5;
        while exp_limit >= 0
            val = integral(f, lo, 10^exp_limit, 'RelTol', 1e-4, 'AbsTol', 1e-6);
            if ~isnan(val), break; end
            exp_limit = exp_limit - 1;
        end
        if isnan(val), val = 1; end
    end
end

% =========================================================================
if strcmp(ch_type, 'Rayl') && strcmp(beamforming_type, 'ANI')
% ── Rayleigh fading, AN beamforming — PDF given by eq. (19) ──────────
    tau  = @(x) 1 + (x .* (1 - alpha)) ./ (alpha .* (k - 1));
    fpdf = @(x) ((tau(x) + SNR*mu_E*(1-alpha)) ./ (SNR*mu_E*alpha)) ...
                .* exp(-x / alpha / SNR / mu_E) .* (tau(x)).^(-k);
    func2 = @(x) myQfunc(func1(x)) .* fpdf(x);
    ail   = safe_integral(func2, x0, inf);

% =========================================================================
elseif strcmp(ch_type, 'Rayl') && strcmp(beamforming_type, 'MRT')
% ── Rayleigh fading, MRT — γ_e ~ Exp(ρβ_e) — PDF given by eq. (23) ──
    fpdf  = @(x) (1 ./ (SNR * mu_E)) .* exp(-x / SNR / mu_E);
    func2 = @(x) myQfunc(func1(x)) .* fpdf(x);
    ail   = safe_integral(func2, x0, inf);

% =========================================================================
elseif strcmp(ch_type, 'Rice') && strcmp(beamforming_type, 'ANI')
% ── Rician fading, AN beamforming — γ_e ≈ Gamma(a_z, b_z) (eq. 34–35)
    a = modelParams(1);   % shape parameter a_z
    b = modelParams(2);   % scale parameter b_z

    y   = sym('y');
    pdf_sym = ((y)/b)^(a-1) .* exp(-(y)/b) / (b * gamma(a));
    fpdf    = matlabFunction(pdf_sym);

    func2 = @(x) myQfunc(func1(x)) .* fpdf(x);
    ail   = safe_integral(func2, x0, inf);

% =========================================================================
elseif strcmp(ch_type, 'Rice') && strcmp(beamforming_type, 'MRT')
% ── Rician fading, MRT — γ_e follows scaled non-central χ² (eq. 39) ──
    noncentrality = modelParams(1);
    sigma_est     = modelParams(2);
    lambda = (noncentrality / sigma_est)^2;
    scale  = sigma_est^2;

    % PDF of non-central chi-square scaled r.v. (r = 2 degrees of freedom)
    r   = 2;
    y   = sym('y');
    pdf_sym    = 0.5 * exp(-(y + lambda)/2) .* besseli(0, sqrt(lambda * y));
    pdf_scaled = (1/scale) * subs(pdf_sym, y, y/scale);
    fpdf       = matlabFunction(pdf_scaled);

    func2 = @(x) (myQfunc(func1(x)) .* fpdf(x));
    ail   = safe_integral(func2, x0, inf);

% =========================================================================
else
    warning('ail_exact: unrecognised channel/beamforming combination. Returning 1.');
    ail = 1;
end

end
