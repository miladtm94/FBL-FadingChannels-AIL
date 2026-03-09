% =========================================================================
% sim_rician_fading.m — AIL Simulation: Rician Fading
% =========================================================================
% Computes the AIL for MRT and AN beamforming over Rician fading channels.
%
% Key difference from Rayleigh: the distribution of Eve's SNR γ_e under
% Rician fading is intractable in closed form (eq. 33). We therefore:
%   1. Generate many γ_e realisations via Monte-Carlo
%   2. Fit a parametric model:
%       - AN  beamforming → Gamma(a_z, b_z)   (eq. 34–35)
%       - MRT beamforming → Rician(s, σ)       (eq. 39)
%   3. Use the fitted model in the analytical AIL expressions
%
% Results are appended to data/processed/FinalResults.mat.
%
% Called by:  simulations/run_all.m
% Depends on: config/sys_params.m, src/*.m
%
% Reference:
%   M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Sec. IV-B.
% =========================================================================

% ── Load system parameters ────────────────────────────────────────────────
sys_params;

% ── Preallocate result arrays ─────────────────────────────────────────────
nN   = length(NVec);
nSNR = length(SNRVec_dB);

deltaBar_Rice_MRT                  = zeros(nN, nSNR);
deltaBar_Rice_MRT_SIM              = zeros(nN, nSNR);
deltaBar_Rice_MRT_APPROX           = zeros(nN, nSNR);
deltaBar_Rice_MRT_APPROX_HighSNR   = zeros(nN, nSNR);

deltaBar_Rice_ANI                  = zeros(nN, nSNR);
deltaBar_Rice_ANI_SIM              = zeros(nN, nSNR);
deltaBar_Rice_ANI_Approx           = zeros(nN, nSNR);
deltaBar_Rice_ANI_Approx_HighSNR   = zeros(nN, nSNR);

% ── Main simulation loop ──────────────────────────────────────────────────
for j = 1:nN
    N = NVec(j);

    for i = 1:nSNR
        SNR = db2pow(SNRVec_dB(i));

        % ── Rician fading: MRT beamforming ───────────────────────────────
        hb = hb_rice;
        [w, U, gammab] = beamforming(hb, 'MRT');
        gammaB = mu_B * SNR * gammab;

        % Fit Rician distribution to γ_e samples
        est_p       = {w, U, 1, mu_E, SNR};
        model_params = dist_estimate(he_rice, 'MRT', est_p);

        % Monte-Carlo simulation
        sim_p = {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR};
        deltaBar_Rice_MRT_SIM(j,i) = ail_simulate(sim_p, 'MRT', he_rice);

        % Exact
        p = {N, M, epsilon, alpha, mu_E, gammaB, SNR, k, Ke_factor, model_params};
        deltaBar_Rice_MRT(j,i) = ail_exact(p, 'Rice', 'MRT');

        % Approximate + High-SNR
        [deltaBar_Rice_MRT_APPROX(j,i), ...
         deltaBar_Rice_MRT_APPROX_HighSNR(j,i)] = ail_approx(p, 'Rice', 'MRT');

        % ── Rician fading: AN beamforming ────────────────────────────────
        [w, U, gammab] = beamforming(hb, 'ANI');
        gammaB = alpha * mu_B * SNR * gammab;

        % Fit Gamma distribution to γ_e samples
        est_p        = {w, U, alpha, mu_E, SNR};
        model_params  = dist_estimate(he_rice, 'ANI', est_p);

        % Monte-Carlo simulation
        sim_p = {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR};
        deltaBar_Rice_ANI_SIM(j,i) = ail_simulate(sim_p, 'ANI', he_rice);

        % Exact
        p = {N, M, epsilon, alpha, mu_E, gammaB, SNR, k, Ke_factor, model_params};
        deltaBar_Rice_ANI(j,i) = ail_exact(p, 'Rice', 'ANI');

        % Approximate + High-SNR
        [deltaBar_Rice_ANI_Approx(j,i), ...
         deltaBar_Rice_ANI_Approx_HighSNR(j,i)] = ail_approx(p, 'Rice', 'ANI');
    end

    fprintf('Rician: blocklength N=%d (%d/%d) done.\n', NVec(j), j, nN);
end

% ── Append to results file ────────────────────────────────────────────────
save(filename, ...
    'deltaBar_Rice_MRT',    'deltaBar_Rice_MRT_SIM', ...
    'deltaBar_Rice_MRT_APPROX', 'deltaBar_Rice_MRT_APPROX_HighSNR', ...
    'deltaBar_Rice_ANI',    'deltaBar_Rice_ANI_SIM', ...
    'deltaBar_Rice_ANI_Approx', 'deltaBar_Rice_ANI_Approx_HighSNR', ...
    '-append');

fprintf('\n[sim_rician_fading] Done. Results appended to %s\n', filename);
