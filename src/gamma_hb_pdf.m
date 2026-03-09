function pdf_val = gamma_hb_pdf(x, params, ch_type)
% GAMMA_HB_PDF  PDF of the main channel power gain γ̃_b.
%
%   pdf_val = GAMMA_HB_PDF(x, params, ch_type)
%
%   Rayleigh fading:
%     γ̃_b ~ Gamma(k, ρβ_b)  →  f(x) = x^{k-1} e^{-x/ρβ_b} / (ρβ_b)^k Γ(k)
%
%   Rician fading:
%     γ̃_b follows a scaled non-central chi-squared distribution with
%     2k degrees of freedom, noncentrality 2kK_b, and scale 2ρβ_b(1+K_b).
%
%   Inputs:
%     x        - evaluation point(s)
%     params   - cell {k, SNR, mu_B, Kb_factor}
%     ch_type  - 'Rayl' or 'Rice'
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Eq. (18), (30).

[k, SNR, mu_B, Kb_factor] = params{:};

if strcmp(ch_type, 'Rayl')
    pdf_val = (x.^(k-1) .* exp(-x ./ (SNR .* mu_B))) ...
              ./ ((SNR .* mu_B).^k .* gamma(k));
else
    scale = 2 * SNR * mu_B * (Kb_factor + 1);
    V     = 2 * k;          % degrees of freedom
    Delta = 2 * k * Kb_factor;  % noncentrality parameter
    % Non-central chi-squared PDF (scaled)
    pdf_val = (1/scale) * 0.5 .* (x/scale ./ Delta).^(V/4 - 1/2) ...
              .* exp(-(x/scale + Delta)/2) ...
              .* besseli(V/2 - 1, sqrt(Delta .* x/scale));
end
end
