function [x,Z,u,v,res,flag0] = hopf_admm_imag_k(A,B,C,mu,la,maxit,tol,k)
%HOPF_ADMM solves
%   min sum_j 0.5*norm(A(:,:,j)*x + (Y+uv')*1i*B(:,j) - C(:,j))^2 + mu*norm(Z,ell_1) + la*norm(v,2), 
%        s.t. Z = Y
%Input:
%   A: space * polynomial basis * time
%   B: space * time
%   C: space * time
%   mu: penalty for sparsity, larger mu implies sparser Z
%   la: penality for norm of v, larger la implies lower-norm v
%   maxit: maximum number of iterations
%   tol: tolerance in stopping rule
%Output:
%   x: coefficients of polynomial basi
%   Z: Sparse coupling matrix (Real-valued)
%   u: Low-rank coupling indicator (Real-valued)
%   v: Low-rank coupling indicator (Real-valued)
%   res: residuals
%   flag0: is convergent
%by JiangnanZhang @Fudan September 20, 2025

% constants
[n,m,t] = size(A);
sqrtnt = sqrt(n*n);
sqrtm = sqrt(m);
sqrtnk = sqrt(n*k);
% initialization
x = zeros(m,1);
Y = zeros(n,n);
Z = zeros(n,n);

u = randn(n,k);
u = u ./ vecnorm(u, 2, 1);
v = randn(n,k);

lambda = zeros(n,n);
rho = min(max(mu*10,1),10000);
tau = 1e-3;
A = reshape(permute(A,[1,3,2]),[n*t,m]);
G = A'*A + eye(m).*tau;
dG = decomposition(G,'chol');
F = real(B)*real(B)' + imag(B)*imag(B)' + eye(n).*tau;
dF0 = decomposition(F + eye(n)*la,'chol');

res = zeros(4,maxit);
flag0 = 1; flag = 1;
% ADMM iteration
for iter = 1:maxit
    x0 = x;
    Z0 = Z;
    u0 = u;
    v0 = v;
    if flag
        dF = decomposition(F + eye(n).*rho,'chol');
        rhosigma = rho*1.618;
        rhotau = 1/(rho+tau);
        murhotau = mu*rhotau;
    end
    % update x
    x = dG\(x.*tau + A'*reshape(C-(Y+u0*v0')*1i*B,[],1));

    % update Y
    W = reshape(A*x,[n,t]);
    Y = (Y.*tau + Z.*rho - lambda + ...
        (real(W-C)-u0*v0'*imag(B))*imag(B)'+...
        + (imag(C-W)-u0*v0'*real(B))*real(B)')/dF;

    % update Z
    Z = (Y.*rho + Z.*tau + lambda).*rhotau;
    Z = sign(Z).*max(abs(Z)-murhotau,0);

    % update rank-1 matrix u*v'
    s0 = (real(W-C)-Y*imag(B))*(imag(B)'*v0)+(imag(C-W)-Y*real(B))*(real(B)'*v0) + tau*u0;
    [u_l,~,v_l] = svd(s0,'econ');
    u = u_l*v_l';
    
    s0 = (real(B)*(imag(C-W)-Y*real(B))'+ imag(B)*(real(W-C)-Y*imag(B))')*u0 + tau*v0;
    v = dF0\s0;

    % update Lambda
    lambda = lambda + (Y-Z).*rhosigma;

    % check convergence
    res(1,iter) = norm(Y-Z,'fro')/sqrtnt;
    res(2,iter) = norm(Z-Z0,'fro')/sqrtnt;
    res(3,iter) = norm(x-x0,2)/sqrtm;
    res(4,iter) = max(norm(u-u0,'fro'),norm(v-v0,'fro'))/sqrtnk;

    if max(res(:,iter)) < tol
        flag0 = 0;
        break
    end

    % adapt penalty
    if res(1,iter) < max(res(2:4,iter))*0.2
        rho = rho*1.1;
        flag = 1;
    elseif res(1,iter) > max(res(2:4,iter))*5
        rho = rho/1.1;
        flag = 1;
    else
        flag = 0;
    end
end
end

