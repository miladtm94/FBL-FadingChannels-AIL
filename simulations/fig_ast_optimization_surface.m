% =========================================================================
% fig_ast_optimization_surface.m — AST Optimisation Surface (Figs. 11–12)
% =========================================================================
% Computes and visualises the Average Secrecy Throughput (AST) as a function
% of power allocation factor α and coding blocklength N for Rayleigh fading.
%
% Identifies the globally optimal (α*, N*) pair by exhaustive grid search,
% then plots the smooth surface and marks the optimum.
%
% Outputs:
%   data/processed/Results_optimalAST.mat
%   results/figures/fig_ast_surface_rayleigh.pdf
%
% Reference:
%   M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Sec. V, Fig. 11.
% =========================================================================
clc; clear; close all;

% ── Paths ─────────────────────────────────────────────────────────────────
repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(genpath(fullfile(repo_root, 'config')));

sys_params;
SNR = db2pow(3);    % fixed transmit SNR = 3 dB

% ── Grid of (α, N) values ────────────────────────────────────────────────
alphaVec = linspace(1e-2, 1,   20);    % 20 power allocation values
NVec_grid= round(linspace(Nmin, Nmax, 100));  % 100 blocklength values

AST_Rayl = zeros(length(alphaVec), length(NVec_grid));

% Pre-compute dispersion at high-SNR asymptote
Vb_inf = channel_dispersion(inf);

% Use a fixed AN beamforming realisation
hb_fixed = normalize(v_B + sigma_B * (randn([1,k]) + 1i*randn([1,k])));
[w, ~, ~] = beamforming(hb_fixed, 'ANI');

fprintf('Computing AST surface (%d × %d grid) ...\n', ...
    length(alphaVec), length(NVec_grid));
tic;

for i = 1:length(alphaVec)
    alpha_i = alphaVec(i);

    for ii = 1:length(NVec_grid)
        N_ii = NVec_grid(ii);
        Rh   = M / N_ii + sqrt(Vb_inf / N_ii) * qfuncinv(epsilon);

        % Find γ_T such that AIL = φ (Rayleigh, AN beamforming)
        syms x_var
        x0_sym   = ((1 + alpha_i * x_var) / (exp(Rh * log(2)))) - 1;
        tau_sym  = 1 + (x0_sym * (1 - alpha_i)) / (alpha_i * (k-1));
        ail_sym  = 1 - exp(-x0_sym / alpha_i / SNR / mu_E) / (tau_sym^(k-1));

        gamma1 = double(vpasolve(ail_sym == 1 - phi, x_var));
        if isempty(gamma1)
            gamma1 = inf;
        end

        % AST: T̄ = (m(1−ε)/N) · Pr{γ̃_b ≥ γ_T}
        % Rayleigh: Pr{γ̃_b ≥ γ_T} = 1 − gammainc(γ_T/ρβ_b, k)  (incomplete Gamma tail)
        AST_Rayl(i, ii) = M * (1-epsilon) * (1 - gammainc(gamma1/SNR/mu_B, k)) / N_ii;
    end

    fprintf('  α = %.3f  (%d/%d)\n', alpha_i, i, length(alphaVec));
end

toc;

% ── Save results ──────────────────────────────────────────────────────────
results_dir = fullfile(repo_root, 'data', 'processed');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save(fullfile(results_dir, 'Results_optimalAST.mat'), 'AST_Rayl', 'alphaVec', 'NVec_grid');

% ── Figures ───────────────────────────────────────────────────────────────
figs_dir = fullfile(repo_root, 'results', 'figures');
if ~exist(figs_dir, 'dir'), mkdir(figs_dir); end

[X, Y] = meshgrid(NVec_grid, alphaVec);

% Raw surface
fig1 = figure('Name', 'Fig 11: AST Surface (raw)');
surf(X, Y, AST_Rayl); colorbar;
xlabel('Blocklength ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Power allocation ($\alpha$)', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('AST $\bar{\mathcal{T}}$ [bpcu]', 'Interpreter', 'latex', 'FontSize', 12);
title('Rayleigh Fading: AST Surface', 'Interpreter', 'latex', 'FontSize', 10);
[row, col]  = find(ismember(AST_Rayl, max(AST_Rayl(:))), 1);
scatter3(X(row,col), Y(row,col), AST_Rayl(row,col), 100, 'red', 'filled', 'o');

% Smoothed surface (for publication)
fig2 = figure('Name', 'Fig 12: AST Surface (smoothed)');
AST_Smooth = smooth_surf(AST_Rayl);
surf(X, Y, AST_Smooth); colorbar;
xlabel('Blocklength ($N$)', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Power allocation ($\alpha$)', 'Interpreter', 'latex', 'FontSize', 12);
zlabel('AST $\bar{\mathcal{T}}$ [bpcu]', 'Interpreter', 'latex', 'FontSize', 12);
title('Rayleigh Fading: AST Surface (smoothed)', 'Interpreter', 'latex', 'FontSize', 10);

[row, col] = find(ismember(AST_Smooth, max(AST_Smooth(:))), 1);
N_opt      = X(row, col);
alpha_opt  = Y(row, col);
scatter3(N_opt, alpha_opt, AST_Smooth(row,col), 100, 'red', 'filled', 'o');

fprintf('\nOptimal solution (smoothed surface):\n');
fprintf('  N*     = %d\n', N_opt);
fprintf('  alpha* = %.4f\n', alpha_opt);
fprintf('  AST*   = %.6f bpcu\n', AST_Smooth(row,col));

saveas(fig2, fullfile(figs_dir, 'fig_ast_surface_rayleigh.pdf'));
fprintf('Saved: results/figures/fig_ast_surface_rayleigh.pdf\n');
