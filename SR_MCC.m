function [x, fail] = SR_MCC(anc, rangems, maxiter, sgm, epsilon)
%Paper: Maximum Correntropy Criterion for Robust TOA-Based Localization in NLOS Environments
%MCC-GTRS Algorithm

%-Inputs
%anc - matrix including sensor positions 
%rangems - measured distance vector
%maxiter - max iteration number
%sgm - initial kernel size
%epsilon - small tolerance for avoiding numerical problems

%-Outputs
%x - location estimate
%fail - true if fail to converge, otherwise false

L = length(anc);
x = zeros(2,1);
p = zeros(L,1);
A1 = [];
b1 = [];
for i = 1:L
    A1 = [A1; -2 * anc(:,i)', 1];
    b1 = [b1; rangems(i)^2 - norm(anc(:,i))^2];
end
D = zeros(3,3);
D(1:2,1:2) = diag([1,1]);
f = [0;0;-0.5];
w = zeros(L,1);

fail = false;

%k: counter for iterations
k = 0;

while 1
    k = k + 1;
    if k > maxiter
        fprintf('have reached the maximum iteration number\n')
        break
    end
    
    x_old = x;
    
    %update p
    for i = 1:L
        p(i) = -exp(-(rangems(i)^2 - norm(x-anc(:,i))^2)^2/(2*sgm^2));
        w(i) = sqrt(-p(i)+epsilon);
    end
    
    %update x, in which GTRS is carried out
    W = diag(w);
    A = W * A1;
    b = W * b1;    
    check1 =  (A'*A)^(1/2) \ D / (A'*A)^(1/2) ;
    
    if sum(sum(isnan(check1))) >= 1 || sum(sum(isinf(check1))) >= 1   
        fail = true;
        return    
    end    
    
    eigen_values = eig( (A'*A)^(1/2) \ D / (A'*A)^(1/2) );
        
    eig_1 = max(eigen_values);

    min_lim = -1/eig_1;

    max_lim = 1e6;

    tol = 1e-3;

    N_iter = 30;

    lambda = bisection_fun(min_lim, max_lim, tol, N_iter, A, D, b, f);
    
    y_hat = (A' * A + lambda * D + 1e-6 * eye(3)) \ (A' * b - lambda * f);

    x = [y_hat(1); y_hat(2)]; % y = [x_1, x_2, ||x||^2, b, b^2]'
    
    if sum(sum(isnan(x))) >= 1 || sum(sum(isinf(x))) >= 1
        
        fail = true;
        return
        
    end
    
    %update sgm
    sgm_old = sgm;
    Error = [];
    for i = 1:L
        Error = [Error; rangems(i)^2 - norm(x - anc(:,i))^2];
    end
    sgm_E = std(Error);
    R = iqr(Error);
    
    sgm = max(1.06*min(sgm_E, R/1.34)*L^(-0.2),sgm_old);
    
    if (norm(x - x_old))<1e-5
        break
    end    
    
end

fprintf('It takes %d iterations due to HQ\n', k)

end

