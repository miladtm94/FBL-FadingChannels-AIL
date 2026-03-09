function Q = generalized_marcum_q(a, b, m)
% GENERALIZED_MARCUM_Q  Compute the generalised Marcum Q-function Q_m(a, b).
%
%   Q = GENERALIZED_MARCUM_Q(a, b, m)
%
%   The generalised Marcum Q-function of order m is defined as (eq. 29):
%
%       Q_m(a, b) = (1/a^{m-1}) ∫_b^∞  y^m · exp(-(y²+a²)/2) · I_{m-1}(ay) dy
%
%   where I_{m-1}(·) is the modified Bessel function of the first kind of
%   order m-1. It is used in the AIL expression for Rician fading (eq. 41)
%   and the non-adaptive AST optimisation for Rician channels.
%
%   Inputs:
%     a  - scalar; noncentrality-related parameter (a ≥ 0)
%     b  - scalar; lower limit of integration (b ≥ 0)
%     m  - scalar; order of the Marcum Q-function (positive integer)
%
%   Output:
%     Q  - scalar; Q_m(a, b) ∈ [0, 1]
%
%   Note:
%     This function uses symbolic integration and may be slow for large m.
%     For m=1, MATLAB's built-in marcumq(a, b) can be used directly.
%
%   Reference:
%     M. Tatar Mamaghani et al., IEEE Trans. Wireless Commun., 2024. Eq. (29).
%     M. K. Simon and M.-S. Alouini, "Digital Communication over Fading
%     Channels," Wiley, 2005.
%
% =========================================================================

y    = sym('y');
func = y^m * exp(-(y^2 + a^2)/2) * besseli(m-1, a*y);
Q    = double((1/a^(m-1)) * int(func, y, b, Inf));

end
