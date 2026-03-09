function Z_smooth = smooth_surf(Z)
% SMOOTH_SURF  2-D Gaussian convolution for surface smoothing.
%
%   Z_smooth = SMOOTH_SURF(Z)
%
%   Applies a 2-D normalised Gaussian kernel to a matrix Z, suitable for
%   removing numerical artefacts from optimisation landscapes (e.g., the
%   AST surface in Figs. 11–12).
%
%   The kernel is an 11×11 Gaussian with standard deviation σ = √5:
%       K(i,j) = exp(-(i²+j²)/10),  normalised to sum to 1.
%
%   Inputs:
%     Z        - m×n matrix; raw surface values (e.g., AST over α × N grid)
%
%   Output:
%     Z_smooth - m×n matrix; smoothed surface (same size as Z)
%
%   Reference:
%     Used in simulations/fig_ast_optimization_surface.m for Fig. 11.
%
% =========================================================================

[xk, yk] = meshgrid(-5:5);
K        = exp(-(xk.^2 + yk.^2) / 10);
K        = K / sum(K(:));                             % normalise kernel

Z_smooth = conv2(Z, K, 'same') ./ conv2(ones(size(Z)), K, 'same');

end
