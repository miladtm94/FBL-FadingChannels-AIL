% =========================================================================
% sim_rayleigh_fading.m — AIL Simulation: Rayleigh Fading
% =========================================================================
% Computes the Average Information Leakage (AIL) for MRT and AN beamforming
% over Rayleigh fading channels across a range of transmit SNRs and
% blocklengths.
%
% Produces (for each combination of N ∈ NVec and SNR ∈ SNRVec_dB):
%   - Exact AIL   (numerical integration of eq. 10)
%   - Approx AIL  (Proposition 1, eq. 20/24)
%   - High-SNR AIL (Sec. IV-C, eq. 44–47)
%   - Monte-Carlo AIL simulation
%
% Results are saved (appended) to data/processed/FinalResults.mat.
%
% Called by:  simulations/run_all.m
% Depends on: config/sys_params.m, src/*.m
%
% Reference:
%   M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Sec. IV-A.
% =========================================================================

% ── Load system parameters ────────────────────────────────────────────────
sys_params;

% ── Preallocate result arrays: [length(NVec) × length(SNRVec_dB)] ────────
nN   = length(NVec);
nSNR = length(SNRVec_dB);

deltaBar_Rayl_MRT                = zeros(nN, nSNR);
deltaBar_Rayl_MRT_SIM            = zeros(nN, nSNR);
deltaBar_Rayl_MRT_Approx         = zeros(nN, nSNR);
deltaBar_Rayl_MRT_Approx_HighSNR = zeros(nN, nSNR);

deltaBar_Rayl_ANI                = zeros(nN, nSNR);
deltaBar_Rayl_ANI_SIM            = zeros(nN, nSNR);
deltaBar_Rayl_ANI_Approx         = zeros(nN, nSNR);
deltaBar_Rayl_ANI_Approx_HighSNR = zeros(nN, nSNR);

% ── Main simulation loop ──────────────────────────────────────────────────
for j = 1:nN
    N = NVec(j);

    for i = 1:nSNR
        SNR = db2pow(SNRVec_dB(i));

        % ── Rayleigh fading: MRT beamforming ─────────────────────────────
        hb = hb_rayl;
        [w, U, gammab] = beamforming(hb, 'MRT');
        gammaB = mu_B * SNR * gammab;

        % Monte-Carlo simulation
        sim_p = {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR};
        deltaBar_Rayl_MRT_SIM(j,i) = ail_simulate(sim_p, 'MRT', he_rayl);

        % Exact (numerical integration)
        p = {N, M, epsilon, alpha, mu_E, gammaB, SNR, k, Ke_factor, 1};
        deltaBar_Rayl_MRT(j,i) = ail_exact(p, 'Rayl', 'MRT');

        % Approximate + High-SNR
        [deltaBar_Rayl_MRT_Approx(j,i), ...
         deltaBar_Rayl_MRT_Approx_HighSNR(j,i)] = ail_approx(p, 'Rayl', 'MRT');

        % ── Rayleigh fading: AN beamforming ──────────────────────────────
        [w, U, gammab] = beamforming(hb, 'ANI');
        gammaB = alpha * mu_B * SNR * gammab;

        % Monte-Carlo simulation
        sim_p = {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR};
        deltaBar_Rayl_ANI_SIM(j,i) = ail_simulate(sim_p, 'ANI', he_rayl);

        % Exact (numerical integration)
        p = {N, M, epsilon, alpha, mu_E, gammaB, SNR, k, Ke_factor, 1};
        deltaBar_Rayl_ANI(j,i) = ail_exact(p, 'Rayl', 'ANI');

        % Approximate + High-SNR
        [deltaBar_Rayl_ANI_Approx(j,i), ...
         deltaBar_Rayl_ANI_Approx_HighSNR(j,i)] = ail_approx(p, 'Rayl', 'ANI');
    end

    fprintf('Rayleigh: blocklength N=%d (%d/%d) done.\n', NVec(j), j, nN);
end

% ── Save results ──────────────────────────────────────────────────────────
save(filename, ...
    'deltaBar_Rayl_MRT',    'deltaBar_Rayl_MRT_SIM', ...
    'deltaBar_Rayl_MRT_Approx', 'deltaBar_Rayl_MRT_Approx_HighSNR', ...
    'deltaBar_Rayl_ANI',    'deltaBar_Rayl_ANI_SIM', ...
    'deltaBar_Rayl_ANI_Approx', 'deltaBar_Rayl_ANI_Approx_HighSNR');

fprintf('\n[sim_rayleigh_fading] Done. Results saved to %s\n', filename);
