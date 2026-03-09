% =========================================================================
% sys_params.m — Centralised System Parameter Configuration
% =========================================================================
% Purpose:
%   Defines all system-level parameters used across simulation scripts.
%   Run this script (via `sys_params`) at the start of any simulation to
%   initialise the shared workspace.
%
% Usage:
%   sys_params          % call from any simulation script
%
% Paper:
%   "Performance Analysis of Finite Blocklength Transmissions Over
%    Wiretap Fading Channels: An Average Information Leakage Perspective"
%   M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024.
% =========================================================================

%% ── Channel / Link Parameters ────────────────────────────────────────────
k       = 4;        % Number of transmit antennas at Alice (MISO system)
mu_B    = 3;        % Large-scale channel gain: Alice → Bob  (β_b = β₀·d_b^{-η})
mu_E    = 1;        % Large-scale channel gain: Alice → Eve  (β_e = β₀·d_e^{-η})

%% ── Rician Fading Parameters ────────────────────────────────────────────
K_factor    = 5;    % Rician K-factor (applies to both links unless overridden)
Kb_factor   = K_factor;
Ke_factor   = K_factor;

zeta_b  = 1;        % Unit-modulus LoS component direction, Bob's link
zeta_e  = 1;        % Unit-modulus LoS component direction, Eve's link

% Derived Rician parameters
sigma_B = sqrt(1 / (2*(Kb_factor + 1)));
v_B     = sqrt(Kb_factor / (Kb_factor + 1)) * zeta_b;  % LoS amplitude, Bob

sigma_E = sqrt(1 / (2*(Ke_factor + 1)));
v_E     = sqrt(Ke_factor / (Ke_factor + 1)) * zeta_e;  % LoS amplitude, Eve

%% ── FBL Coding Parameters ────────────────────────────────────────────────
epsilon = 1e-3;     % Target decoding error probability at Bob (ε)
M       = 100;      % Number of confidential information bits per packet (m)
N       = 400;      % Default coding blocklength in channel uses
Nmin    = 1;        % Minimum allowed blocklength (optimisation lower bound)
Nmax    = 1e3;      % Maximum allowed blocklength (optimisation upper bound)

%% ── Beamforming / Power Allocation ──────────────────────────────────────
alpha   = 0.7;      % Power fraction allocated to information signal (0 < α ≤ 1)
                    % Remaining (1−α) is distributed to Artificial Noise (AN)
                    % Special case: α = 1 → MRT beamforming (no AN)

%% ── Optimisation ─────────────────────────────────────────────────────────
phi     = 1e-4;     % Maximum tolerated AIL constraint (φ) for AST optimisation

%% ── Transmit SNR ─────────────────────────────────────────────────────────
SNR     = 1;        % Default transmit SNR ρ = P/σ² (linear scale; see also db2pow)

%% ── Simulation Settings ──────────────────────────────────────────────────
simNum  = 1e4;      % Number of Monte-Carlo channel realisations
L       = 1e3;      % Number of channel realisations for adaptive scheme

%% ── Pre-generate Channel Realisations ───────────────────────────────────
% Eve's channel realisations (used by Rayleigh and Rician simulations)
he_rayl = zeros(simNum, k);
he_rice = zeros(simNum, k);

for i = 1:simNum
    htilde       = (randn([1, k]) + 1i * randn([1, k]));
    he_rayl(i,:) = sqrt(0.5) * htilde;                  % CN(0, I_k), normalised
    he_rice(i,:) = v_E + sigma_E * htilde;               % Rician
end

% Row-wise normalisation (unit norm per realisation)
he_rayl = normalize(he_rayl, 2);
he_rice = normalize(he_rice, 2);

% Bob's channel (a single fixed realisation used per simulation run)
hb_rayl = normalize(sqrt(0.5) * (randn([1, k]) + 1i * randn([1, k])));

hb_rice = v_B + sigma_B * (randn([1, k]) + 1i * randn([1, k]));
hb_rice = normalize(hb_rice);
