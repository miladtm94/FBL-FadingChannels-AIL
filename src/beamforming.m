function [w, U, gammaB] = beamforming(hb, type)
% BEAMFORMING  Generate beamforming vectors for AN or MRT strategy.
%
%   [w, U, gammaB] = BEAMFORMING(hb, type)
%
%   Inputs:
%     hb    - 1×k complex row vector; normalised small-scale fading coefficient
%             of the Alice-Bob (main) channel.
%     type  - string; beamforming strategy:
%               'ANI' — Artificial Noise (AN) beamforming (see Sec. IV-A)
%               'MRT' — Maximum Ratio Transmission (MRT), i.e., α = 1
%
%   Outputs:
%     w      - k×1 complex unit vector; information beamforming vector
%              (principal eigenvector of hb'*hb, points towards Bob)
%     U      - k×(k-1) complex matrix; AN precoding matrix
%              (orthonormal basis for the null-space of hb)
%              For MRT, U is a zero matrix (no AN).
%     gammaB - scalar; Bob's normalised channel power gain |hb·w|²
%
%   System model (eq. 14 in the paper):
%     x = sqrt(α)·w·s + sqrt((1-α)/(k-1))·U·v
%       where s is the unit-power information symbol and v ~ CN(0, I_{k-1})
%
%   Reference:
%     M. Tatar Mamaghani et al., "Performance Analysis of Finite Blocklength
%     Transmissions Over Wiretap Fading Channels," IEEE Trans. Wireless
%     Commun., 2024. Eq. (14)–(17).
%
% =========================================================================

k = length(hb);

if strcmp(type, 'ANI')
    % ── AN Beamforming (Sec. IV-A) ─────────────────────────────────────
    % Eigen-decompose H = hb'*hb to separate the signal and null spaces.
    H = hb' * hb;                          % k×k rank-1 Hermitian matrix
    [V, D] = eig(H);
    [~, ind] = sort(diag(D), 'descend');   % sort eigenvalues descending
    Vs = V(:, ind);

    w = Vs(:, 1);                          % principal eigenvector → Bob
    U = Vs(:, 2:end);                      % remaining k-1 eigenvectors → null space

elseif strcmp(type, 'MRT')
    % ── MRT Beamforming ───────────────────────────────────────────────
    % Maximise SNR at Bob; no AN is used (α = 1 implicitly).
    w = hb' ./ norm(hb);                   % matched filter / MRT vector
    U = zeros(k, k-1);                     % no AN precoder

else
    error('beamforming: unknown type "%s". Use ''ANI'' or ''MRT''.', type);
end

% ── Bob's received SNR (normalised, without large-scale gain or transmit SNR)
gammaB  = hb * (w * w') * hb';
gammaB  = (gammaB + gammaB') / 2;         % enforce real-valued output

end
