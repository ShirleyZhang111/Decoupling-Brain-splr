function [coefficient, sparse, u, v, v_hat, res, mse_real,cos_theta] = generate_coupling_components(inCfg)
% Individual nonlinear model
% Input: 
%   TC : space * time
%   T : time points for fitting
%   coupling mechanism : bipartite (Bipartite Conjugate Coupling) /
%   amplitude (Amplitude-based Coupling)
% Output:
%   coefficient: coefficients of polynomial basis
%   Sparse: sparse component
%   u,v: u*v' is the rank-k component
%   v_hat: embedded v
%   res: residuals
%   mse_real : mean squared error of data fitting
%   cos_theta : cosine similarity between the subspaces spanned by U from different runs
%   rank: rank of the low-rank component 

% input configuration
data_path = inCfg.data_path;
file_name = inCfg.file_name;
save_path = inCfg.save_path;
if isfield(inCfg,'rank');rank = inCfg.rank;else;rank = 0;end
if isfield(inCfg,'maxit');maxit = inCfg.maxit;else;maxit = 500;end
if isfield(inCfg,'tol');tol = inCfg.tol;else;tol = 1e-4;end
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end
if isfield(inCfg,'xi');xi = inCfg.xi;else;xi = single(0.1*length(T));end
if isfield(inCfg,'mu');mu = inCfg.mu;else;mu = single(0.1*length(T));end
if isfield(inCfg,'mechanism');mechanism = inCfg.mechanism;else;mechanism = 'bipartite';end 
if ~contains(file_name,'.mat');file_name = [file_name '.mat'];end

% start processing
load([data_path file_name],'TC');% N*T
TC = TC(:,6:end-5);
if isempty(T);T = 5:size(TC,2)-4;end

hilbert_data = double(hilbert(zscore(TC'))).';clearvars TC
B = hilbert_data(:,T);
t0 = round(size(B,2)/5*4);
C = conv2(hilbert_data(:,[T(1:4)-4 T T(end-3:end)+4]),fliplr([1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280]),'valid');%.*-1i
% C = conv2(hilbert_data(:,[T(1)-1 T T(end)+1]),[1/2 0 -1/2],'valid');
hilbert_data = permute(hilbert_data(:,T),[1,3,2]);
A = [ones(size(hilbert_data)),...
    hilbert_data,conj(hilbert_data),...
    hilbert_data.*hilbert_data,hilbert_data.*conj(hilbert_data),conj(hilbert_data.*hilbert_data),...
    hilbert_data.*hilbert_data.*hilbert_data,hilbert_data.*hilbert_data.*conj(hilbert_data),hilbert_data.*conj(hilbert_data.*hilbert_data),conj(hilbert_data.*hilbert_data.*hilbert_data)];
clearvars hilbert_data

mse_real = []; cos_theta = [];
if strcmp(mechanism, 'bipartite') == 1
    if length(mu)==1 % calculate results with optimal sparsity parameter
        [coefficient, sparse, u, v, res, isconverge] = hopf_admm_imag_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
        [~, S0, V0] = svd(B(:,1:t0)'); s = diag(S0); v_hat = (V0'*v).*s;

        if isconverge;n_iter = length(find(res(1,:)>0));else;n_iter = maxit;end
        disp(['Subject ',file_name(1:end-4),' finished, use time ',num2str(n_iter*toc/60),' mins, residual: ' num2str(max(res(:,n_iter)))]);
        % Save the result to a file
        save([save_path file_name(1:end-4) '_fitting_' mechanism '.mat'], 'coefficient', 'sparse','u','v','res','isconverge','mu','xi','tol','rank');
        
    else % find optimal sparsity parameter
        for i_mu = length(mu):-1:1
            for i_xi = length(xi):-1:1
                [coefficient, sparse, u, v, res, ~] = hopf_admm_imag_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
                U0(:,:,1) = u;
                delta_Z = C(:,t0+1:end);
                coupling_mat = (sparse + u*v').*1i;
                dZ = construct_data(B(:,t0+1:end),coefficient,coupling_mat);
                mse_real(1,i_mu,i_xi) = norm(real(delta_Z-dZ),'fro');
                for iter = 2:10
                    [coefficient, sparse, u, v, res, ~] = hopf_admm_imag_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
                    U0(:,:,iter) = u; 
                end   
                for iter1 = 1:9
                    for iter2 = iter1+1:10
                        Cos_theta0(iter1,iter2) = norm(U0(:,:,iter1)'*U0(:,:,iter2),2);
                    end
                end
                cos_theta(1,i_mu,i_xi) = mean(Cos_theta0(Cos_theta0 ~= 0));
            end
        end   
    end
elseif strcmp(mechanism, 'amplitude') == 1
    if length(mu)==1 % calculate results with optimal sparsity parameter
        [coefficient, sparse, u, v, res, isconverge] = hopf_admm_real_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
        [~, S0, V0] = svd(B(:,1:t0)'); s = diag(S0); v_hat = (V0'*v).*s;
        if isconverge;n_iter = length(find(res(1,:)>0));else;n_iter = maxit;end
        disp(['Subject ',file_name(1:end-4),' finished, use time ',num2str(n_iter*toc/60),' mins, residual: ' num2str(max(res(:,n_iter)))]);
        % Save the result to a file
        save([save_path file_name(1:end-4) '_fitting_' mechanism '.mat'], 'coefficient', 'sparse','u','v', 'res','isconverge','mu','xi','tol','rank');
        
    else % find optimal sparsity parameter
        for i_mu = length(mu):-1:1
            for i_xi = length(xi):-1:1
                [coefficient, sparse, u, v, res, ~] = hopf_admm_real_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
                U0(:,:,1) = u;
                delta_Z = C(:,t0+1:end);
                coupling_mat = sparse + u*v';
                dZ = construct_data(B(:,t0+1:end),coefficient,coupling_mat);
                mse_real(1,i_mu,i_xi) = norm(real(delta_Z-dZ),'fro');
                for iter = 2:10
                    [coefficient, sparse, u, v, res, ~] = hopf_admm_real_k(A(:,:,1:t0),B(:,1:t0),C(:,1:t0),mu,xi,maxit,tol,rank);
                    U0(:,:,iter) = u; 
                end   
                for iter1 = 1:9
                    for iter2 = iter1+1:10
                        Cos_theta0(iter1,iter2) = norm(U0(:,:,iter1)'*U0(:,:,iter2),2);
                    end
                end
                cos_theta(1,i_mu,i_xi) = mean(Cos_theta0(Cos_theta0 ~= 0));
            end
        end   
    end
end
