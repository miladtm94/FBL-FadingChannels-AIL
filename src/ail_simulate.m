function ail_mc = ail_simulate(sim_params, beamforming_type, He)
% AIL_SIMULATE  Monte-Carlo Average Information Leakage estimation.
%
%   ail_mc = AIL_SIMULATE(sim_params, beamforming_type, He)
%
%   Estimates the AIL by averaging the per-realisation information leakage
%   δ(γ_e) over many independent channel realisations of Eve's channel.
%   This function provides the ground-truth Monte-Carlo reference used to
%   validate the analytical expressions in the paper.
%
%   Per-realisation leakage (from eq. 9):
%       δ(γ_e) = Q(√(N/V_e) · [C_s(γ_b, γ_e) − R₀])
%
%   Inputs:
%     sim_params        - cell array:
%                           {gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR}
%                         gammaB — Bob's received SNR (fixed per simulation)
%                         w      — k×1 info beamforming vector
%                         U      — k×(k-1) AN precoding matrix
%     beamforming_type  - string; 'ANI' or 'MRT'
%     He                - simNum×k complex matrix; Eve's channel realisations
%                         (each row is one normalised realisation h_e)
%
%   Output:
%     ail_mc            - scalar; empirical average δ̄ = mean(δ(γ_e))
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Eq. (9).
%
% =========================================================================

[simNum, k] = size(He);
[gammaB, N, epsilon, M, alpha, w, U, mu_E, SNR] = sim_params{:};

Rs  = M / N;
Vb  = channel_dispersion(gammaB);
R   = Rs + sqrt(Vb ./ N) .* qfuncinv(epsilon);  % effective codeword rate R₀

delta = zeros(1, simNum);

for n = 1:simNum
    he = He(n, :);

    % ── Compute Eve's instantaneous SNR γ_e ─────────────────────────
    if strcmp(beamforming_type, 'ANI')
        % AN beamforming SINR at Eve (eq. 17)
        gammaE = alpha * he * w * ...
                 (((1-alpha)/(k-1)) * he * (U*U') * he' + 1/mu_E/SNR)^(-1) ...
                 * w' * he';
    else
        % MRT: Eve sees only direct leakage (eq. 21)
        gammaE = mu_E * SNR * he * (w * w') * he';
    end
    gammaE = (gammaE + gammaE') / 2;   % enforce real-valued

    % ── Per-realisation leakage δ(γ_e) = Q(⋅) ───────────────────────
    Cs       = log2((1 + gammaB) ./ (1 + gammaE));      % IBL secrecy capacity
    Ve       = channel_dispersion(gammaE);
    delta(n) = qfunc(sqrt(N ./ Ve) .* (Cs - R));
end

ail_mc = mean(delta);

end
