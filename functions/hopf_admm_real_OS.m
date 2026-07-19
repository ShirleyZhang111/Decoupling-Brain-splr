function [x,Z,res,flag,flag0] = hopf_admm_real_OS(A,B,C,mu,maxit,tol)
%HOPF_ADMM solves (OS: only-sparse model, which means ran == 1)
%   min sum_j 0.5*norm(A(:,:,j)*x + Z*B(:,j) - C(:,j))^2 + mu*norm(Z,ell_1)
%Input:
%   A: space * polynomial basis * time
%   B: space * time
%   C: space * time
%   mu: penalty for sparsity, larger mu implies sparser Z
%   maxit: maximum number of iterations
%   tol: tolerance in stopping rule
%Output:
%   x: coefficients of polynomial basis
%   Z: coupling matrix (Real-valued)
%   res: residuals
%by JiangnanZhang @Fudan September 20, 2025

% constants
[n,m,t] = size(A);
sqrtnt = sqrt(n*t);
sqrtm = sqrt(m);
% initialization
x = zeros(m,1);
Y = zeros(n,n);
Z = zeros(n,n);
lambda = zeros(n,n);
rho = min(max(mu*10,1),10000);
tau = 1e-3;
A = reshape(permute(A,[1,3,2]),[n*t,m]);
G = A'*A + eye(m).*tau;
dG = decomposition(G,'chol');
F = real(B)*real(B)' + imag(B)*imag(B)' + eye(n).*tau;
res = zeros(2,maxit);
flag0 = 1; flag = 1;
% ADMM iteration
for iter = 1:maxit
    x0 = x;
    Z0 = Z;
    if flag
        dF = decomposition(F + eye(n).*rho,'chol');
        rhosigma = rho*1.618;
        rhotau = 1/(rho+tau);
        murhotau = mu*rhotau;
    end
    % update x
    x = dG\(x.*tau + A'*reshape(C-Y*B,[],1));
    
    % update Y
    W = reshape(A*x,[n,t]);
    Y = (Y.*tau + Z.*rho - lambda + real(W-C)*real(B)' + imag(W-C)*imag(B)')/dF;
    % Y = (Y.*tau + Z.*rho - lambda + (C-reshape(A*x,[n,t]))*B')/dF;

    % update Z
    Z = (Y.*rho + Z.*tau + lambda).*rhotau;
    Z = sign(Z).*max(abs(Z)-murhotau,0);

    % update Lambda
    lambda = lambda + (Y-Z).*rhosigma;

    % check convergence
    res(1,iter) = norm(Y-Z,'fro')/sqrtnt;
    res(2,iter) = norm(x-x0,2)/sqrtm;
    res(3,iter) = norm(Z-Z0,'fro')/sqrtnt;

    if max(res(:,iter)) < tol
        flag0 = 0;
        break
    end
    
    % adapt penalty
    if res(1,iter) < mean(res(2:3,iter))*0.1
        rho = rho*0.618;
        flag = 1;
    elseif res(1,iter) > mean(res(2:3,iter))*10
        rho = rho/0.618;
        flag = 1;
    else
        flag = 0;
    end
end

end

