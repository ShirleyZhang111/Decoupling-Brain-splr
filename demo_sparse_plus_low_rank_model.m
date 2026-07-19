% Sparse-plus-low-rank model
% This code implements the sparse-plus-low-rank model for fMRI data analysis

% Load rs-fMRI signals (available datasets: HCP, HCPex, HCP-voxel, UKB)
data_path = './data/';
file_name = 'HCPex_rest_data.mat';
save_path = './results/';

% Set configuration parameters for nonlinear model
cfg.rank = 1; % rank of low-rank component rank = 0,1,2,3,...
cfg.maxit = 2000;
cfg.tol = 1e-12;

% coupling mechanism:
% bipartite: Bipartite conjugate coupling mechanism 
% amplitude: Amplitude-based coupling mechanism
cfg.mechanism = 'bipartite'; 

cfg.T = 101:600; 
cfg.mu = 0.1*length(cfg.T);
cfg.xi = 100;
cfg.data_path = data_path;
cfg.file_name = file_name;
cfg.save_path = save_path;

% Estimate Model Parameters from rs-fMRI data using nonlinear model
[coeff, sparse, u, v, ~, ~, ~, ~] = generate_coupling_components(cfg);

% Nonzeros pattren of sparse coupling matrix
figure; spy(sparse);

% Calculate the correlation between u-factor and gradient (individual, cerebral cortex)
load('data\HCP360_gradient.mat')
u = u.*sign(u(1,:));
[r1,p1] = corr(u(1:360),gradient_HCP360);

% Calculate the correlation between u-factor and GSR_R2 (individual)
R2_GSR2 = calc_GSR_R2(cfg);
[r2,p2] = corr(u,R2_GSR2);