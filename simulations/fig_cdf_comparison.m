% =========================================================================
% fig_cdf_comparison.m — CDF Validation of γ_e Gamma Approximation (Fig. 2)
% =========================================================================
% Validates the Gamma distribution approximation of Eve's SNR γ_e under
% AN beamforming over Rician fading channels.
%
% Compares F_{γ_e}(x) (Monte-Carlo CDF) with F_Z(x) (fitted Gamma CDF)
% for three transmit SNR values: ρ ∈ {-3, 0, 3} dB.
%
% Reference: M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Fig. 2.
% =========================================================================
clc; clear; close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'config')));

sys_params;
simNum_cdf = 1e5;   % large sample count for smooth CDF estimate

% Generate Rician channel realisations
rng('default');
He_cdf  = v_E + sigma_E * (randn([simNum_cdf,k]) + 1i*randn([simNum_cdf,k]));
hb_cdf  = v_B + sigma_B * (randn([1,k]) + 1i*randn([1,k]));
[w, U, ~] = beamforming(normalize(hb_cdf), 'ANI');

SNR_vec_cdf = db2pow([-3, 0, 3]);
xVec        = linspace(0, 7, 1e3);

% Compute CDFs
cdf_sim    = zeros(length(xVec), length(SNR_vec_cdf));
cdf_approx = zeros(length(xVec), length(SNR_vec_cdf));

for j = 1:length(SNR_vec_cdf)
    SNR_j = SNR_vec_cdf(j);

    % Collect γ_e samples
    gammaE_samples = zeros(simNum_cdf, 1);
    for n = 1:simNum_cdf
        he = He_cdf(n,:);
        gE = alpha * he * w * (((1-alpha)/(k-1)) * he*(U*U')*he' + 1/mu_E/SNR_j)^(-1) * w' * he';
        gammaE_samples(n) = (gE + conj(gE))/2;
    end

    % Fit Gamma distribution
    pd = fitdist(gammaE_samples, 'Gamma');
    Xest = makedist('Gamma', 'a', pd.a, 'b', pd.b);

    for ii = 1:length(xVec)
        cdf_sim(ii, j)    = mean(gammaE_samples < xVec(ii));
        cdf_approx(ii, j) = cdf(Xest, xVec(ii));
    end
end

% ── Plot ──────────────────────────────────────────────────────────────────
curveColours = {"#0072BD", "#A2142F", "#77AC30"};
fig = figure('Name', 'Fig. 2: CDF Validation');

for j = 1:length(SNR_vec_cdf)
    plot(xVec, cdf_sim(:,j), '-', 'Color', curveColours{j}, 'LineWidth', 2); hold on;
    plot(xVec, cdf_approx(:,j), 's', 'Color', curveColours{j}, 'LineWidth', 1.5, ...
        'MarkerIndices', 1:round(length(xVec)/30):length(xVec), 'MarkerSize', 6);
end

xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 12);
title('$F_{\gamma_e}(x)$ [Sim.] vs. $F_Z(x)$ [Approx.]', ...
    'Interpreter', 'latex', 'FontSize', 10);
legend('Sim.', 'Approx.', 'FontSize', 9, 'Location', 'southeast', ...
    'Interpreter', 'latex', 'FontName', 'Times New Roman');
grid on;

figs_dir = fullfile(repo_root, 'results', 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end
saveas(fig, fullfile(figs_dir, 'fig_cdf_comparison.pdf'));
fprintf('Saved: results/figures/fig_cdf_comparison.pdf\n');
