% =========================================================================
% run_all.m — Master Simulation Runner
% =========================================================================
% Runs all simulation experiments and reproduces Figs. 3–4 from the paper:
%   "Performance Analysis of Finite Blocklength Transmissions Over
%    Wiretap Fading Channels: An Average Information Leakage Perspective"
%
% Pipeline:
%   1. Simulate Rayleigh fading (sim_rayleigh_fading.m) → saves to .mat
%   2. Simulate Rician fading  (sim_rician_fading.m)   → appends to .mat
%   3. Load results and generate publication-quality figures
%
% Estimated runtime: 2–8 hours (depends on hardware; simNum = 1e4)
%
% Outputs:
%   data/processed/FinalResults.mat   — all AIL result arrays
%   results/figures/fig_rayleigh_ail.pdf
%   results/figures/fig_rician_ail.pdf
%
% Usage:
%   >> run_all
%
% =========================================================================
clc; clear; close all;

% ── Paths ─────────────────────────────────────────────────────────────────
repo_root  = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'config')));
cd(fileparts(mfilename('fullpath')));   % cd to simulations/

% ── Sweep parameters ──────────────────────────────────────────────────────
SNRVec_dB = linspace(-20, 10, 30);     % transmit SNR range [dB]
NVec      = [5e2, 5e3, 5e4];           % blocklengths to evaluate

% Results file
results_dir  = fullfile(repo_root, 'data', 'processed');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
filename = fullfile(results_dir, 'FinalResults.mat');

% ── Step 1: Rayleigh fading simulation ───────────────────────────────────
fprintf('\n=== [1/2] Rayleigh Fading Simulation ===\n');
sim_rayleigh_fading;

clearvars -except filename SNRVec_dB NVec repo_root;

% ── Step 2: Rician fading simulation ─────────────────────────────────────
fprintf('\n=== [2/2] Rician Fading Simulation ===\n');
sim_rician_fading;

% ── Step 3: Load results and plot ─────────────────────────────────────────
fprintf('\n=== Generating Figures ===\n');
sys_params;   % reload system parameters for plot labels

markerIdx      = 1:floor(length(SNRVec_dB)/20):length(SNRVec_dB);
curveColours   = {"#0072BD", "#A2142F", "#77AC30", '#7E2F8E', '#000000'};
curveMarkers   = {'<', 'p', 's', 'o'};

load(filename);

figs_dir = fullfile(repo_root, 'results', 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

% ── Figure 1: Rayleigh fading ─────────────────────────────────────────────
fig1 = figure('Name', 'Fig 3: AIL vs SNR (Rayleigh)');

layout = tiledlayout(2, 1);

ax1 = nexttile;
for j = 1:length(NVec)
    semilogy(ax1, SNRVec_dB, deltaBar_Rayl_ANI(j,:), ...
        'Color', curveColours{1}, 'LineStyle', 'none', ...
        'Marker', curveMarkers{1}, 'LineWidth', 2, 'MarkerIndices', markerIdx);
    hold on;
    semilogy(ax1, SNRVec_dB, deltaBar_Rayl_ANI_Approx(j,:), ...
        'Color', curveColours{2}, 'LineStyle', '-', 'LineWidth', 2);
    semilogy(ax1, SNRVec_dB, deltaBar_Rayl_ANI_Approx_HighSNR(j,:), ...
        'Color', curveColours{3}, 'LineStyle', ':', 'LineWidth', 2);
end
set(ax1, 'YScale', 'log');
legend('Exact', 'Approx.', 'High-SNR', 'Interpreter', 'latex', ...
    'FontSize', 8, 'Location', 'best', 'FontName', 'Times New Roman');
title(sprintf('AN beamforming ($\\alpha$ = %0.1f)', alpha), ...
    'FontSize', 10, 'Interpreter', 'latex');
grid on; grid minor;

ax2 = nexttile;
for j = 1:length(NVec)
    semilogy(ax2, SNRVec_dB, deltaBar_Rayl_MRT(j,:), ...
        'Color', curveColours{1}, 'LineStyle', 'none', ...
        'Marker', curveMarkers{2}, 'LineWidth', 2, 'MarkerIndices', markerIdx);
    hold on;
    semilogy(ax2, SNRVec_dB, deltaBar_Rayl_MRT_Approx(j,:), ...
        'Color', curveColours{2}, 'LineStyle', '-', 'LineWidth', 2);
    semilogy(ax2, SNRVec_dB, deltaBar_Rayl_MRT_Approx_HighSNR(j,:), ...
        'Color', curveColours{3}, 'LineStyle', ':', 'LineWidth', 2);
end
set(ax2, 'YScale', 'log');
legend('Exact', 'Approx.', 'High-SNR', 'Interpreter', 'latex', ...
    'FontSize', 8, 'Location', 'best', 'FontName', 'Times New Roman');
title('MRT beamforming', 'FontSize', 10, 'Interpreter', 'latex');
grid on; grid minor;

linkaxes([ax1, ax2], 'x');
xticklabels(ax1, {});
layout.TileSpacing = 'compact';
ylabel(layout, 'Average Information Leakage ($\bar{\delta}$)', ...
    'FontSize', 12, 'Interpreter', 'latex');
xlabel(layout, 'Transmit SNR ($\rho$) [dB]', 'FontSize', 12, 'Interpreter', 'latex');

saveas(fig1, fullfile(figs_dir, 'fig_rayleigh_ail.pdf'));
fprintf('Saved: results/figures/fig_rayleigh_ail.pdf\n');

% ── Figure 2: Rician fading ───────────────────────────────────────────────
fig2 = figure('Name', 'Fig 4: AIL vs SNR (Rician)');

layout2 = tiledlayout(2, 1);

ax3 = nexttile;
for j = 1:length(NVec)
    semilogy(ax3, SNRVec_dB, deltaBar_Rice_ANI(j,:), ...
        'Color', curveColours{1}, 'LineStyle', 'none', ...
        'Marker', curveMarkers{1}, 'LineWidth', 2, 'MarkerIndices', markerIdx);
    hold on;
    semilogy(ax3, SNRVec_dB, deltaBar_Rice_ANI_Approx(j,:), ...
        'Color', curveColours{2}, 'LineStyle', '-', 'LineWidth', 2);
    semilogy(ax3, SNRVec_dB, deltaBar_Rice_ANI_Approx_HighSNR(j,:), ...
        'Color', curveColours{3}, 'LineStyle', ':', 'LineWidth', 2);
end
set(ax3, 'YScale', 'log');
legend('Exact', 'Approx.', 'High-SNR', 'Interpreter', 'latex', ...
    'FontSize', 8, 'Location', 'southwest', 'FontName', 'Times New Roman');
title(sprintf('AN beamforming ($\\alpha$ = %0.1f)', alpha), ...
    'FontSize', 10, 'Interpreter', 'latex');
grid on; grid minor;

ax4 = nexttile;
for j = 1:length(NVec)
    semilogy(ax4, SNRVec_dB, deltaBar_Rice_MRT(j,:), ...
        'Color', curveColours{1}, 'LineStyle', 'none', ...
        'Marker', curveMarkers{2}, 'LineWidth', 2, 'MarkerIndices', markerIdx);
    hold on;
    semilogy(ax4, SNRVec_dB, deltaBar_Rice_MRT_APPROX(j,:), ...
        'Color', curveColours{2}, 'LineStyle', '-', 'LineWidth', 2);
    semilogy(ax4, SNRVec_dB, deltaBar_Rice_MRT_APPROX_HighSNR(j,:), ...
        'Color', curveColours{3}, 'LineStyle', ':', 'LineWidth', 2);
end
set(ax4, 'YScale', 'log');
legend('Exact', 'Approx.', 'High-SNR', 'Interpreter', 'latex', ...
    'FontSize', 8, 'Location', 'northeast', 'FontName', 'Times New Roman');
title('MRT beamforming', 'FontSize', 10, 'Interpreter', 'latex');
grid on; grid minor;

linkaxes([ax3, ax4], 'x');
xticklabels(ax3, {});
layout2.TileSpacing = 'compact';
ylabel(layout2, 'Average Information Leakage ($\bar{\delta}$)', ...
    'FontSize', 12, 'Interpreter', 'latex');
xlabel(layout2, 'Transmit SNR ($\rho$) [dB]', 'FontSize', 12, 'Interpreter', 'latex');

saveas(fig2, fullfile(figs_dir, 'fig_rician_ail.pdf'));
fprintf('Saved: results/figures/fig_rician_ail.pdf\n');

fprintf('\n=== All simulations complete. ===\n');
