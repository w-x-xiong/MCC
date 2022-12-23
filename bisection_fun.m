function [lambda] = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, g)

lmbd_lower = min_lim;
lmbd_upper = max_lim;
lmbd_mid = (lmbd_lower+lmbd_upper)/2;

Func_Y = @(lmbd) ( ((A)'*(A)+lmbd*D) + 1e-6 * eye(3) ) \ ( (A')*b - lmbd*g );
Func_Phi = @(lmbd) ((Func_Y(lmbd))')*D*Func_Y(lmbd) + 2*(g')*Func_Y(lmbd);

maxiter = 0;

while abs(Func_Phi(lmbd_mid)) > tol
    if (Func_Phi(lmbd_mid)*Func_Phi(lmbd_upper)) < 0
        lmbd_lower = lmbd_mid;
    else
        lmbd_upper = lmbd_mid;
    end
    lmbd_mid = (lmbd_lower+lmbd_upper)/2;
    maxiter = maxiter + 1;
    if maxiter > N_iter
        break
    end
end
lambda = lmbd_mid;

end