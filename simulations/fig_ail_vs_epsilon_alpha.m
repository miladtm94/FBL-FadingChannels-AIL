% =========================================================================
% fig_ail_vs_epsilon_alpha.m — AIL vs Error Probability & Alpha (Figs. 5–6)
% =========================================================================
% Plots the AIL as a function of the decoding error probability ε for
% different blocklengths N and power allocation values α, for both
% Rayleigh and Rician fading channels.
%
% Shows that relaxing the reliability constraint (larger ε) reduces AIL,
% revealing the fundamental reliability–secrecy trade-off in FBL systems.
%
% Reference: M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Figs. 5–6.
% =========================================================================
clc; clear; close all;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'config')));

sys_params;

epsVec   = logspace(-10, -1, 20);  % error probability sweep
NVec     = [3e2, 5e2, 1e3];         % blocklength values
alphaVec = [0.3, 0.7, 1.0];         % power split values

nE  = length(epsVec);
nN  = length(NVec);
nA  = length(alphaVec);

deltaBar_Rayl_Exact       = zeros(nE, nN, nA);
deltaBar_Rayl_Approx      = zeros(nE, nN, nA);
deltaBar_Rice_Exact       = zeros(nE, nN, nA);
deltaBar_Rice_Approx      = zeros(nE, nN, nA);

for ii = 1:nN
    N_ii = NVec(ii);
    for iii = 1:nE
        eps_iii = epsVec(iii);
        for i = 1:nA
            alpha_i = alphaVec(i);

            if alpha_i ~= 1
                % ── AN beamforming ────────────────────────────────────
                [w, U, gammab] = beamforming(hb_rayl, 'ANI');
                gammaB = alpha_i * mu_B * SNR * gammab;
                p = {N_ii, M, eps_iii, alpha_i, mu_E, gammaB, SNR, k, Ke_factor, 1};
                [deltaBar_Rayl_Approx(iii,ii,i), ~] = ail_approx(p, 'Rayl', 'ANI');
                deltaBar_Rayl_Exact(iii,ii,i)       = ail_exact(p,  'Rayl', 'ANI');

                [w, U, gammab] = beamforming(hb_rice, 'ANI');
                gammaB = alpha_i * mu_B * SNR * gammab;
                est_p  = {w, U, alpha_i, mu_E, SNR};
                mp     = dist_estimate(he_rice, 'ANI', est_p);
                p = {N_ii, M, eps_iii, alpha_i, mu_E, gammaB, SNR, k, Ke_factor, mp};
                [deltaBar_Rice_Approx(iii,ii,i), ~] = ail_approx(p, 'Rice', 'ANI');
                deltaBar_Rice_Exact(iii,ii,i)       = ail_exact(p,  'Rice', 'ANI');
            else
                % ── MRT beamforming ───────────────────────────────────
                [w, U, gammab] = beamforming(hb_rayl, 'MRT');
                gammaB = mu_B * SNR * gammab;
                p = {N_ii, M, eps_iii, alpha_i, mu_E, gammaB, SNR, k, Ke_factor, 1};
                [deltaBar_Rayl_Approx(iii,ii,i), ~] = ail_approx(p, 'Rayl', 'MRT');
                deltaBar_Rayl_Exact(iii,ii,i)       = ail_exact(p,  'Rayl', 'MRT');

                [w, U, gammab] = beamforming(hb_rice, 'MRT');
                gammaB = mu_B * SNR * gammab;
                est_p  = {w, U, 1, mu_E, SNR};
                mp     = dist_estimate(he_rice, 'MRT', est_p);
                p = {N_ii, M, eps_iii, alpha_i, mu_E, gammaB, SNR, k, Ke_factor, mp};
                [deltaBar_Rice_Approx(iii,ii,i), ~] = ail_approx(p, 'Rice', 'MRT');
                deltaBar_Rice_Exact(iii,ii,i)       = ail_exact(p,  'Rice', 'MRT');
            end
        end
    end
    fprintf('N = %d (%d/%d) done.\n', NVec(ii), ii, nN);
end

% ── Figures ───────────────────────────────────────────────────────────────
curveColours  = {"#0072BD", "#A2142F", "#77AC30"};
curveMarkers  = {'<', 'p', 's'};
curveLines    = {'-', '--', '-.'};
markerIdx     = 0:2:length(epsVec); markerIdx(1) = 1;

figs_dir = fullfile(repo_root, 'results', 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

for fading = {'Rayleigh', 'Rician'}
    fig = figure('Name', sprintf('AIL vs epsilon (%s)', fading{1}));
    layout = tiledlayout(1, 3);

    for j = 1:nA
        ax = nexttile;
        for i = 1:nN
            if strcmp(fading{1}, 'Rayleigh')
                exact_data  = deltaBar_Rayl_Exact(:,i,j);
                approx_data = deltaBar_Rayl_Approx(:,i,j);
            else
                exact_data  = deltaBar_Rice_Exact(:,i,j);
                approx_data = deltaBar_Rice_Approx(:,i,j);
            end
            loglog(ax, epsVec, exact_data,  'Color', curveColours{j}, ...
                'LineStyle', curveLines{j}, 'LineWidth', 2);
            hold on;
            loglog(ax, epsVec, approx_data, 'Color', curveColours{j}, ...
                'Marker', curveMarkers{j}, 'LineStyle', 'none', ...
                'LineWidth', 1.5, 'MarkerIndices', markerIdx);
        end
        title(ax, sprintf('$\\alpha = %.2f$', alphaVec(j)), ...
            'Interpreter', 'latex', 'FontSize', 12);
        legend(ax, 'Exact', 'Approx.', 'Interpreter', 'latex', ...
            'FontSize', 8, 'Location', 'northeast');
        grid(ax, 'on'); grid(ax, 'minor');
    end

    linkaxes(findall(fig, 'Type', 'axes'), 'x');
    ylabel(layout, 'Average Information Leakage ($\bar{\delta}$)', ...
        'Interpreter', 'latex', 'FontSize', 12);
    xlabel(layout, 'Decoding Error Probability ($\varepsilon$)', ...
        'Interpreter', 'latex', 'FontSize', 12);
    title(layout, sprintf('%s Fading', fading{1}), ...
        'Interpreter', 'latex', 'FontSize', 10);

    fname = sprintf('fig_ail_vs_epsilon_%s.pdf', lower(fading{1}));
    saveas(fig, fullfile(figs_dir, fname));
    fprintf('Saved: results/figures/%s\n', fname);
end
