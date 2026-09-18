function [x, y, lambda, ok] = wsrls_gtrs(A1, b1, D, f, weights)
% WSRLS_GTRS  Globally solve the weighted squared-range least-squares
% problem
%
%   min_y  || diag(sqrt(weights)) (A1*y - b1) ||_2^2
%   s.t.   y' * D * y + 2 * f' * y = 0,           y = [x; alpha] in R^{H+1}
%
% via the generalized trust region subproblem (GTRS) machinery of
% More (1993): the global solution is y(lambda*) = (C + lambda* D)^{-1}
% (c - lambda* f), where lambda* is the unique root, on
% I = (-1/chi1, inf), chi1 = largest generalized eigenvalue of the
% matrix pencil (D, C), of the strictly-decreasing secular function
% psi(lambda) = y(lambda)' D y(lambda) + 2 f' y(lambda).
%
% This function factors out and generalizes (to arbitrary dimension H)
% the GTRS-solving logic that, in the original prototype, was written
% inline inside SR_MCC.m for H = 2 only. It is used both for the
% unweighted warm start and for every half-quadratic iteration of
% SR_MCC.m.
%
% Corrections relative to the original inline code:
%   (i)   chi1 is obtained directly as the largest generalized
%         eigenvalue of the pencil (D, C) via eig(D, C), instead of
%         forming the explicit matrix square root C^{1/2} and computing
%         eig(C^{-1/2} D C^{-1/2}); this avoids one matrix square root
%         and one redundant recomputation of the same product, and is
%         numerically more stable when C is ill-conditioned;
%   (ii)  a relative (trace-scaled) Tikhonov regularizer is used in
%         place of the fixed absolute constant 1e-6, so behaviour does
%         not depend on the (arbitrary) physical units of the problem;
%   (iii) the function reports failure (ok = false) rather than
%         propagating NaN/Inf silently.

L = size(A1, 1);
n = size(D, 1);          % n = H + 1
H = n - 1;

if any(~isfinite(weights)) || sum(weights) < 1e-13 * numel(weights)
    x = nan(H, 1); y = nan(n, 1); lambda = NaN; ok = false;
    return
end

w = sqrt(max(weights, 0));
A = A1 .* w(:);
b = b1 .* w(:);

C = A' * A;
c = A' * b;

scale = max(trace(C) / n, 1e-300);
reg = 1e-10 * scale;

gev = eig(D, C + reg * eye(n));   % generalized eigenvalues of (D, C)
gev = sort(real(gev));
chi1 = gev(end);

if chi1 > 0
    min_lim = -1 / chi1 * (1 - 1e-9);
else
    min_lim = -1e-8 * max(scale, 1);
end
max_lim = max(1, abs(min_lim)) * 10;   % verified / expanded inside bisection_fun

tol_bisect = 1e-10 * max(1, scale);
N_iter = 100;

[lambda, ok] = bisection_fun(min_lim, max_lim, tol_bisect, N_iter, A, D, b, f);
if ~ok
    x = nan(H, 1); y = nan(n, 1);
    return
end

M = C + lambda * D + reg * eye(n);
y = M \ (c - lambda * f);
if any(~isfinite(y))
    x = nan(H, 1); ok = false;
    return
end
x = y(1:H);

end
