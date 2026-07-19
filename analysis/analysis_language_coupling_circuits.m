% Hemispheric Specialization of Sparse Coupling Circuits during Language Processing
% Group-level sparse-plus-rank-1 model
load('../results/group_language_Coup_rank1.mat')

S = mean(abs(Sparse_Coup_language), 3);
S = S - diag(diag(S));

% Right hemisphere: nodes 1–180; 
S_right = S(1:180, 1:180);

% Set a threshold and retain only strong connections
threshold = prctile(abs(S_right(:)), 95);  % Retain the top 5% of connections
A = abs(S_right) > threshold;              % Binary adjacency matrix
W = S_right .* A;                          % Weighted connectivity matrix
L = 5;       % Maximum cycle length (recommended range: 4–8)
K = 50;      % Number of strongest cycles to return
tol = 0;     % Use 1e-12 to remove numerically negligible edges

top_length7_right = top_cycles_exactL_sumabs(W, L, K, tol);

for i = 1:numel(top_length7_right)
    fprintf('#%d len=%d score=mean|w|=%g nodes: ', ...
        i, top_length7_right(i).len, top_length7_right(i).score / L);
    fprintf('%d ', top_length7_right(i).nodes);
    fprintf('\n  edge |w|: ');
    fprintf('%g ', top_length7_right(i).edgeAbs);
    fprintf('\n');
end

load('../data/HCPex_LabelID.mat');
% Sort labels according to their region indices
LabelID = sortrows(LabelID, 2);
% Store the anatomical labels of nodes in each right-hemisphere cycle
path_right = cell(K, L+2);

for i = 1:K
    path_right(i, 1:L) = LabelID(top_length7_right(i).nodes, 3)';
    path_right(i, L + 2) = num2cell(top_length7_right(i).score / L);
end
path_right(:, L + 1) = path_right(:, 1);

% Left Hemisphere: nodes 181–360
S_left = S(181:360,181:360);
threshold = prctile(abs(S_left(:)), 95);  
A = abs(S_left) > threshold;  
W = S_left .* A; 
top_length7_left = top_cycles_exactL_sumabs(W, L, K, tol);
for i = 1:numel(top_length7_left )
    fprintf('#%d len=%d score=mean|w|=%g nodes: ', i, top_length7_left(i).len, top_length7_left(i).score/L);
    fprintf('%d ', top_length7_left(i).nodes);
    fprintf('\n  edge |w|: ');
    fprintf('%g ', top_length7_left(i).edgeAbs);
    fprintf('\n');
end

path_left = cell(K, L+2); 
for i = 1:K
    path_left(i,1:L) = LabelID(top_length7_left(i).nodes + 180,3)';
    path_left(i,L+2) = num2cell(top_length7_left(i).score/L);
end
path_left(:,L+1) = path_left(:,1);
%% Calculate Jaccard index
for k = 1:K
    set_L(k,:) = top_length7_left(k).nodes;
    set_R(k,:) = top_length7_right(k).nodes;
end

intersection = intersect(reshape(set_L,1,[]), reshape(set_R,1,[]));
unions = union(reshape(set_L,1,[]), reshape(set_R,1,[]));
% Jaccard index
J = length(intersection) / length(unions);