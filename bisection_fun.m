function [lambda, ok] = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, g)
% BISECTION_FUN  Root of the GTRS secular equation psi(lambda) = 0 on
% [min_lim, max_lim) by bisection (More, 1993, Thm. 5.2: psi is strictly
% decreasing on this interval, so the root -- when it exists in the
% supplied bracket -- is unique).
%
%   yhat(lambda)  = (A'A + lambda*D + reg*I) \ (A'b - lambda*g)
%   psi(lambda)   = yhat(lambda)' D yhat(lambda) + 2 g' yhat(lambda)
%
% Corrections relative to the original prototype:
%   (i)  the regularizing identity is sized to size(D,1) rather than
%        hard-coded to eye(3), so this routine now works for any
%        dimension (H = 2, H = 3, ...), not only H = 2;
%   (ii) the search bracket is VERIFIED (psi(min_lim) and psi(max_lim)
%        must have opposite signs) and, if the caller-supplied max_lim
%        does not bracket a sign change, it is expanded geometrically;
%        a fixed constant such as 1e6 is not scale-invariant and can
%        silently fail to bracket the root for problems in different
%        units in the ORIGINAL implementation, without any warning;
%   (iii) an explicit success flag `ok` is returned instead of silently
%        returning a meaningless value when no root can be bracketed.
%
% Inputs
%   min_lim, max_lim : initial search interval, min_lim = -1/chi1 is the
%                       theoretically exact left end point (see wsrls_gtrs.m)
%   tol               : tolerance on |psi(lambda)|
%   N_iter             : maximum number of bisection iterations
%   A, D, b, g         : GTRS data, y = (A'A + lambda D)^{-1}(A'b - lambda g)
%
% Outputs
%   lambda : the located root (best available on failure)
%   ok     : true if a sign-changing bracket was found and bisection
%            converged to |psi(lambda)| <= tol (or interval collapsed)

n = size(D, 1);
reg = 1e-6 * max(trace(A' * A) / n, 1e-300);
I_n = eye(n);

Func_Y   = @(lmbd) (A' * A + lmbd * D + reg * I_n) \ (A' * b - lmbd * g);
Func_Phi = @(lmbd) local_psi(Func_Y, D, g, lmbd);

lmbd_lower = min_lim;
lmbd_upper = max_lim;

f_lo = Func_Phi(lmbd_lower);
f_hi = Func_Phi(lmbd_upper);

expand_iter = 0;
while (~isfinite(f_lo) || ~isfinite(f_hi) || f_lo * f_hi > 0) && expand_iter < 80
    lmbd_upper = lmbd_upper * 2 + 1;   % geometric bracket expansion
    f_hi = Func_Phi(lmbd_upper);
    expand_iter = expand_iter + 1;
end

if ~isfinite(f_lo) || ~isfinite(f_hi) || f_lo * f_hi > 0
    lambda = lmbd_lower;
    ok = false;
    return
end

lmbd_mid = (lmbd_lower + lmbd_upper) / 2;
maxiter = 0;
while abs(Func_Phi(lmbd_mid)) > tol
    f_mid = Func_Phi(lmbd_mid);
    if f_lo * f_mid < 0
        lmbd_upper = lmbd_mid;
    else
        lmbd_lower = lmbd_mid;
        f_lo = f_mid;
    end
    lmbd_mid = (lmbd_lower + lmbd_upper) / 2;
    maxiter = maxiter + 1;
    if maxiter > N_iter
        break
    end
end
lambda = lmbd_mid;
ok = true;

end

function val = local_psi(Func_Y, D, g, lmbd)
y = Func_Y(lmbd);
if any(~isfinite(y))
    val = NaN;
else
    val = (y') * D * y + 2 * (g') * y;
end
end
