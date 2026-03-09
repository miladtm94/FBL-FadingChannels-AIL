function [ail, ail_highsnr] = ail_approx(params, ch_type, beamforming_type)
% AIL_APPROX  Approximate Average Information Leakage (Proposition 1).
%
%   [ail, ail_highsnr] = AIL_APPROX(params, ch_type, beamforming_type)
%
%   Computes the closed-form Laplace/saddle-point approximation of the AIL
%   (Proposition 1, eq. 11) and its high-SNR counterpart (Sec. IV-C) for
%   four supported channel / beamforming combinations.
%
%   Inputs:
%     params            - cell array:
%                           {N, M, epsilon, alpha, mu_E, gammaB, SNR, k,
%                            Ke_factor, modelParams}
%                         where modelParams is:
%                           Rayleigh : scalar 1 (unused)
%                           Rician AN: [shape, scale] (Gamma fit of γ_e)
%                           Rician MRT: [noncentrality, sigma] (Rician fit)
%     ch_type           - string; 'Rayl' or 'Rice'
%     beamforming_type  - string; 'ANI' (AN beamforming) or 'MRT'
%
%   Outputs:
%     ail        - scalar; approximate AIL δ̄ at the given SNR
%     ail_highsnr- scalar; high-SNR asymptote of AIL
%
%   Closed-form expressions (paper eqs.):
%     Rayleigh AN  (eq. 20):  δ̄ = exp(−x₀/α/ρ/β_e) / τ(x₀)^{k−1}
%     Rayleigh MRT (eq. 24):  δ̄ = exp(−x₀/ρ/β_e)
%     Rician AN    (eq. 36):  δ̄ = 1 − F_Gamma(x₀; a_z, b_z)
%     Rician MRT   (eq. 41):  δ̄ = Q₁(√(2kK_e), √(2(1+K_e)x₀/ρ/β_e))
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024.
%
% =========================================================================

[N, M, epsilon, alpha, mu_E, gammaB, SNR, k, Ke_factor, modelParams] = params{:};

% Evaluation points under normal and high-SNR modes
x0         = eval_point(gammaB, 'normal', {N, epsilon, M});
x0_highsnr = eval_point(gammaB, 'high',   {N, epsilon, M});

% ── Helper: τ(x) for Rayleigh AN beamforming (eq. 19) ──────────────────
tau     = @(x) 1 + (x * (1 - alpha)) / (alpha * (k - 1));

% =========================================================================
if strcmp(ch_type, 'Rayl') && strcmp(beamforming_type, 'ANI')
% ── Case: Rayleigh fading, AN beamforming (Sec. IV-A, Case I) ────────
% Approximate AIL: δ̄ ≈ exp(−x₀/α/ρβ_e) / τ(x₀)^{k−1}  (eq. 20)

    ail        = exp(-x0         / alpha / SNR / mu_E) / (tau(x0))^(k-1);
    ail_highsnr= exp(-x0_highsnr / alpha / SNR / mu_E) / (tau(x0_highsnr))^(k-1);

% =========================================================================
elseif strcmp(ch_type, 'Rayl') && strcmp(beamforming_type, 'MRT')
% ── Case: Rayleigh fading, MRT beamforming (Sec. IV-A, Case II) ──────
% Approximate AIL: δ̄ ≈ exp(−x₀/ρβ_e)  (eq. 24)

    ail        = exp(-x0         / (mu_E * SNR));
    ail_highsnr= exp(-x0_highsnr / (mu_E * SNR));

% =========================================================================
elseif strcmp(ch_type, 'Rice') && strcmp(beamforming_type, 'ANI')
% ── Case: Rician fading, AN beamforming (Sec. IV-B, Case I) ──────────
% γ_e is approximated by a Gamma distribution Z ~ G(a_z, b_z) (eq. 35)
% Approximate AIL: δ̄ ≈ 1 − F_Z(x₀)  (eq. 36)

    shape_est = modelParams(1);     % a_z (shape)
    scale_est = modelParams(2);     % b_z (scale)
    pd = makedist('Gamma', 'a', shape_est, 'b', scale_est);

    ail        = 1 - cdf(pd, x0);
    ail_highsnr= 1 - cdf(pd, x0_highsnr);

% =========================================================================
elseif strcmp(ch_type, 'Rice') && strcmp(beamforming_type, 'MRT')
% ── Case: Rician fading, MRT beamforming (Sec. IV-B, Case II) ────────
% γ_e follows a scaled non-central chi-square distribution (eq. 39–40)
% Approximate upper-bound AIL: δ̄ ≈ Q₁(√(2kK_e), √(2(1+K_e)x₀/ρβ_e)) (eq. 41)

    noncentrality = modelParams(1);
    sigma_est     = modelParams(2);
    lambda = (noncentrality / sigma_est)^2;
    scale  = sigma_est^2;

    % Use Symbolic Math for the noncentral chi-square PDF (r=2 degrees of freedom)
    y   = sym('y');
    r   = 2;
    pdf_sym      = 0.5 * exp(-(y + lambda)/2) .* ((y/lambda)^(r/4 - 1/2)) .* ...
                   besseli(r/2 - 1, sqrt(lambda * y));
    pdf_scaled   = (1/scale) * subs(pdf_sym, y, y/scale);
    pdf_func     = matlabFunction(pdf_scaled);

    ail        = 1 - integral(pdf_func, 0, x0,         'RelTol', 1e-4, 'AbsTol', 1e-6);
    ail_highsnr= 1 - integral(pdf_func, 0, x0_highsnr, 'RelTol', 1e-4, 'AbsTol', 1e-6);

% =========================================================================
else
    warning('ail_approx: unrecognised ch_type/beamforming combination. Returning 1.');
    ail        = 1;
    ail_highsnr= 1;
end

end
