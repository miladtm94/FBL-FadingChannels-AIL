% =========================================================================
% quick_start.m — Minimal Reproducible Demo
% =========================================================================
% Demonstrates how to compute the Average Information Leakage (AIL) for a
% single parameter configuration, without running the full heavy simulations.
%
% This script:
%   1. Sets up a MISO wiretap channel (k=4 antennas, Rayleigh fading)
%   2. Generates AN beamforming vectors
%   3. Computes and compares Exact, Approximate, and Monte-Carlo AIL
%   4. Prints a results table and shows a simple SNR-sweep plot
%
% Expected runtime: < 2 minutes
%
% Usage:
%   >> cd examples
%   >> quick_start
% =========================================================================
clc; clear; close all;

% ── Add source paths ──────────────────────────────────────────────────────
repo_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(repo_root, '..', 'src')));
addpath(genpath(fullfile(repo_root, '..', 'config')));

%% ── System Parameters (override sys_params defaults as needed) ───────────
k         = 4;         % Alice's antennas
epsilon   = 1e-3;      % Target decoding error probability
M         = 100;       % Information bits per packet
N         = 500;       % Coding blocklength
alpha     = 0.7;       % AN power fraction (0 < α ≤ 1)
mu_B      = 3;         % Alice-Bob path gain β_b
mu_E      = 1;         % Alice-Eve path gain β_e
simNum    = 2000;      % Monte-Carlo realisations (reduced for speed)

% SNR sweep
SNR_dB_vec = linspace(-15, 10, 20);

%% ── Generate Channel Realisations ────────────────────────────────────────
% Bob's channel (single fixed realisation)
hb = normalize(sqrt(0.5) * (randn([1,k]) + 1i*randn([1,k])));

% Eve's channel ensemble (simNum realisations)
He_rayl = normalize(sqrt(0.5) * (randn([simNum,k]) + 1i*randn([simNum,k])), 2);

%% ── Beamforming ──────────────────────────────────────────────────────────
[w, U, gammab] = beamforming(hb, 'ANI');

%% ── AIL Computation over SNR range ──────────────────────────────────────
ail_exact_vec   = zeros(1, length(SNR_dB_vec));
ail_approx_vec  = zeros(1, length(SNR_dB_vec));
ail_mc_vec      = zeros(1, length(SNR_dB_vec));

fprintf('%-10s %-14s %-14s %-14s\n', 'SNR [dB]', 'Exact', 'Approx.', 'Monte-Carlo');
fprintf('%s\n', repmat('-', 1, 55));

for i = 1:length(SNR_dB_vec)
    SNR    = db2pow(SNR_dB_vec(i));
    gammaB = alpha * mu_B * SNR * gammab;

    % Parameter cell array
    p = {N, M, epsilon, alpha, mu_E, gammaB, SNR, k, 5, 1};

    % Exact AIL (numerical integration)
    ail_exact_vec(i)  = ail_exact(p, 'Rayl', 'ANI');

    % Approximate AIL (Proposition 1, closed form)
    [ail_approx_vec(i), ~] = ail_approx(p, 'Rayl', 'ANI');

    % Monte-Carlo simulation
    sim_p = {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR};
    ail_mc_vec(i) = ail_simulate(sim_p, 'ANI', He_rayl);

    fprintf('%-10.1f %-14.4e %-14.4e %-14.4e\n', ...
        SNR_dB_vec(i), ail_exact_vec(i), ail_approx_vec(i), ail_mc_vec(i));
end

%% ── Plot ─────────────────────────────────────────────────────────────────
figure('Name', 'Quick Start: AIL vs SNR (Rayleigh AN beamforming)');

semilogy(SNR_dB_vec, ail_mc_vec,    'ko',  'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Monte-Carlo (sim.)');
hold on;
semilogy(SNR_dB_vec, ail_exact_vec, 'b--', 'LineWidth', 2, ...
    'DisplayName', 'Exact (num. integration)');
semilogy(SNR_dB_vec, ail_approx_vec,'r-',  'LineWidth', 2, ...
    'DisplayName', 'Approx. (Proposition 1)');
grid on; grid minor;

xlabel('Transmit SNR $\rho$ [dB]', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Average Information Leakage $\bar{\delta}$', 'Interpreter', 'latex', 'FontSize', 12);
title(sprintf('AIL vs SNR: AN beamforming, Rayleigh fading\n$k=%d$, $N=%d$, $m=%d$, $\\epsilon=10^{%d}$, $\\alpha=%.1f$', ...
    k, N, M, log10(epsilon), alpha), 'Interpreter', 'latex', 'FontSize', 10);

lgd = legend('Interpreter', 'latex', 'FontSize', 10, 'Location', 'southwest');
lgd.FontName = 'Times New Roman';

fprintf('\nPlot displayed. Simulation complete.\n');
