function x0 = eval_point(gammaB, mode, params)
% EVAL_POINT  Compute the Laplace approximation evaluation point x₀.
%
%   x₀ = EVAL_POINT(gammaB, mode, params)
%
%   In Proposition 1 (paper eq. 11–12), the AIL approximation is
%
%       δ̄ ≈ 1 − F_{γ_e}(x₀)
%
%   where the evaluation point x₀ is defined as:
%
%       x₀ = max{ (1 + γ_b) / 2^{R₀} − 1, 0 }
%
%   and R₀ = √(V_b/N)·Q⁻¹(ε) + m/N is the effective rate overhead.
%
%   Two modes are supported:
%     'normal'  — uses the true dispersion V(γ_b) (finite SNR)
%     'high'    — uses the asymptotic dispersion V(∞) = log₂²(e)
%                  (high-SNR approximation, eq. 42–43)
%
%   Inputs:
%     gammaB  - scalar; Bob's received SNR γ_b (linear)
%     mode    - string; 'normal' or 'high'
%     params  - cell array {N, epsilon, M}:
%                 N       — coding blocklength
%                 epsilon — target decoding error probability
%                 M       — number of information bits
%
%   Output:
%     x₀      - scalar ≥ 0; evaluation point for the AIL approximation
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024.
%     Proposition 1 (eq. 12) and Sec. IV-C (eq. 42–43).
%
% =========================================================================

[N, epsilon, M] = params{:};
Rs = M / N;                         % secrecy rate (bits per channel use)

if strcmp(mode, 'high')
    % ── High-SNR asymptote: V(∞) = log₂²(e) ─────────────────────────
    Vb   = channel_dispersion(inf);
    Rinf = Rs + sqrt(Vb ./ N) .* qfuncinv(epsilon);
    x0   = max(gammaB / (2^Rinf), 0);

else
    % ── Normal (finite-SNR) mode ─────────────────────────────────────
    Vb = channel_dispersion(gammaB);
    R  = Rs + sqrt(Vb ./ N) .* qfuncinv(epsilon);
    x0 = max((1 + gammaB) / (2^R) - 1, 0);
end

end
