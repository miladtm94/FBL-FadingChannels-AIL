function V = channel_dispersion(gamma)
% CHANNEL_DISPERSION  Compute the channel dispersion V(γ).
%
%   V = CHANNEL_DISPERSION(gamma)
%
%   The channel dispersion characterises the stochastic variation of the
%   channel and appears in the normal approximation of the FBL coding rate.
%   For an AWGN channel with SNR γ it is defined as (eq. 8 in the paper):
%
%       V(γ) = log₂²(e) · γ(γ+2) / (γ+1)²
%
%   This function is vectorised: gamma may be a scalar, vector, or matrix.
%
%   Inputs:
%     gamma - SNR value(s) γ (linear, not dB).  Pass Inf for the high-SNR
%             asymptote V(∞) = log₂²(e).
%
%   Output:
%     V     - Channel dispersion value(s), same size as gamma.
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Eq. (8).
%     Y. Polyanskiy et al., IEEE Trans. Inf. Theory, 2010.
%
% =========================================================================

V = log2(exp(1)).^2 .* (1 - (1 + gamma).^(-2));

end
