function model_params = dist_estimate(He_rice, beamforming_type, est_params)
% DIST_ESTIMATE  Estimate the parametric distribution of Eve's SNR γ_e
%                under Rician fading via moment matching / MLE.
%
%   model_params = DIST_ESTIMATE(He_rice, beamforming_type, est_params)
%
%   For Rician fading channels, the exact distribution of Eve's SNR γ_e
%   under AN beamforming is intractable (eq. 33). This function fits a
%   parametric model by generating many realisations of γ_e and fitting:
%     - AN beamforming  → Gamma distribution G(a_z, b_z)  [eq. 34–35]
%     - MRT beamforming → Rician distribution Rice(s, σ)
%
%   The Gamma approximation accuracy is demonstrated in Fig. 2 of the paper.
%
%   Inputs:
%     He_rice           - simNum×k complex matrix; Rician Eve channel realizations
%     beamforming_type  - string; 'ANI' or 'MRT'
%     est_params        - cell array {w, U, alpha, mu_E, SNR}
%                           w     — k×1 info beamforming vector
%                           U     — k×(k-1) AN precoding matrix
%                           alpha — power allocation factor
%                           mu_E  — Eve's large-scale gain β_e
%                           SNR   — transmit SNR ρ
%
%   Output:
%     model_params      - 1×2 vector:
%                           [a_z, b_z]    for AN  (Gamma shape, scale)
%                           [s, sigma]    for MRT (Rician params)
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024.
%     Sec. IV-B, eq. (34)–(35) and surrounding discussion.
%
% =========================================================================

[simNum, k] = size(He_rice);
[w, U, alpha, mu_E, SNR] = est_params{:};

gammaE_samples = zeros(simNum, 1);

for n = 1:simNum
    he = He_rice(n, :);

    if strcmp(beamforming_type, 'ANI')
        % Eve's SINR under AN beamforming (eq. 17)
        gammaE = alpha * he * w * ...
                 (((1-alpha)/(k-1)) * he * (U*U') * he' + 1/mu_E/SNR)^(-1) ...
                 * w' * he';
    else
        % Eve's SNR under MRT (eq. 21)
        gammaE = mu_E * SNR * he * (w * w') * he';
    end

    gammaE_samples(n) = (gammaE + conj(gammaE)) / 2;  % real-valued
end

% ── Fit parametric distribution ──────────────────────────────────────────
if strcmp(beamforming_type, 'ANI')
    % Fit Gamma distribution to γ_e samples (eq. 34)
    pd = fitdist(gammaE_samples, 'Gamma');
    model_params = [pd.a, pd.b];        % [shape, scale]

else
    % Fit Rician distribution to √(γ_e) samples
    pd = fitdist(sqrt(gammaE_samples), 'Rician');
    model_params = [pd.s, pd.sigma];    % [noncentrality, sigma]
end

end
