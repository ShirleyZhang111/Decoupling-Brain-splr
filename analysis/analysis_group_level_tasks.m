% Group-level analysis of task-state coupling component (S+rank-1)
% 7 tasks
load('../results/group_emotion_Coup_rank1.mat')
load('../results/group_gambling_Coup_rank1.mat')
load('../results/group_language_Coup_rank1.mat')
load('../results/group_motor_Coup_rank1.mat')
load('../results/group_relational_Coup_rank1.mat')
load('../results/group_social_Coup_rank1.mat')
load('../results/group_wm_Coup_rank1.mat')

%% Analysis of Generative rank-1 component U
% U_wm,U_relational,U_language,U_motor,U_social,U_emotion,U_gambling
% each task-state U: regions*subjects (379*395)
% U_rest: regions*rank*subjects (379*1*395)
% task and rest: same 395 subjects
[regions,rank,subjects] = size(U_wm);

U_all = [reshape(U_wm,regions,[]),...
    reshape(U_motor,regions,[]),...
    reshape(U_relational,regions,[]),...
    reshape(U_language,regions,[]),...
    reshape(U_social,regions,[]),...
    reshape(U_emotion,regions,[]),...
    reshape(U_gambling,regions,[])];

% random 400 samples: U_7tasks 379*400
rng('shuffle');               
idx = randperm(size(U_all,2), 400); 
U_7tasks = U_all(:,idx);

% t-SNE visualization for task/rest u-factor
load('../results/group_rest_Coup_rank1.mat')
U_rest = reshape(U_rest,379,[]);
Result = [U_7tasks,U_rest];
data = abs(Result)';
data(isnan(data)) = 0;
Y_tsne = tsne(data, 'NumDimensions', 3, 'Perplexity', 50, 'Verbose', 1);
figure;
hold on;
scatter(Y_tsne(1:400, 1), Y_tsne(1:400, 2),20, 'filled', 'MarkerFaceColor', [37 125 139]/255);
scatter(Y_tsne(401:end, 1), Y_tsne(401:end, 2),20, 'filled','MarkerFaceColor',[189  119 149]/255);
axis square;

%  1000-fold cross-validation SVM performance
group1 = U_7tasks';
group2 = U_rest';
X = [group1; group2];
X = abs(X);
X(isnan(X)) = 0;
Y = [ones(size(group1, 1), 1); 
     2*ones(size(group2, 1), 1)];
accuracy = zeros(10,1);
for fold = 1:10
    rng('shuffle');    
    cv = cvpartition(Y, 'HoldOut', 0.2);
    idx_train = training(cv);
    idx_test = test(cv);
    X_train = X(idx_train, :);
    Y_train = Y(idx_train);
    X_test = X(idx_test, :);
    Y_test = Y(idx_test);
    svm_model = fitcecoc(X_train, Y_train, ...
        'Learners', 'svm', ...
        'Coding', 'onevsone', ...  
        'Verbose', 1);
   
    Y_pred = predict(svm_model, X_test);
    accuracy(fold,1) = sum(Y_pred == Y_test) / length(Y_test);
end
%% Analysis of Generative sparse component S (7 tasks)
[m,n,num_samples] = size(Sparse_Coup_wm);
C1 = cat(3,Sparse_Coup_wm,Sparse_Coup_relational,Sparse_Coup_language,Sparse_Coup_motor,...
    Sparse_Coup_social,Sparse_Coup_emotion,Sparse_Coup_gambling);
X_flat = zeros(num_samples, n*m);
for i = 1:num_samples*7
    X_flat(i, :) = reshape(C1(:,:,i)-diag(diag(C1(:,:,i))), 1, n*m);
end
% dim: PCA dimensionality
dim = 800;
[coeff_7tasks, X, latent_7tasks, ~, explained_7tasks] = pca(X_flat, 'NumComponents', dim);
X(isnan(X)) = 0;
Y = kron((1:7)',ones(num_samples,1));

%  1000-fold cross-validation SVM performance
acc_7tasks = zeros(1000,1);
for fold = 1:1000
    cv = cvpartition(Y, 'HoldOut', 0.2);
    idx_train = training(cv);
    idx_test = test(cv);
    
    X_train = X(idx_train, :);
    Y_train = Y(idx_train);
    X_test = X(idx_test, :);
    Y_test = Y(idx_test);
    
    svm_model = fitcecoc(X_train, Y_train, ...
        'Learners', 'svm', ...
        'Coding', 'onevsall', ...  
        'Verbose', 1);       
    Y_pred = predict(svm_model, X_test);    
    accuracy = sum(Y_pred == Y_test) / length(Y_test);
    acc_7tasks(fold,1) = accuracy;
end

%% Analysis of Generative sparse component S (sub tasks)
% eg.working-memory state 0bk/2bk
load('../results/group_wm_subtasks_Coup_rank1.mat');

[m,n,num_samples] = size(Sparse_Coup_0bk);
C1 = cat(3,Sparse_Coup_0bk,Sparse_Coup_2bk);
X_flat = zeros(num_samples*2, n*m);
for i = 1:num_samples*2
    X_flat(i, :) = reshape(C1(:,:,i)-diag(diag(C1(:,:,i))), 1, n*m);
end
dim = num_samples*2-1;
[coeff, X, latent, ~, explained] = pca(X_flat, 'NumComponents', dim);
X(isnan(X)) = 0;
Y = kron((1:2)',ones(num_samples,1));

%  1000-fold cross-validation SVM performance
folds = 1000;
acc_wm_subtasks = zeros(folds,1);
for fold = 1:folds
    cv = cvpartition(Y, 'HoldOut', 0.2);
    idx_train = training(cv);
    idx_test = test(cv);
    
    X_train = X(idx_train, :);
    Y_train = Y(idx_train);
    X_test = X(idx_test, :);
    Y_test = Y(idx_test);
    
    svm_model = fitcecoc(X_train, Y_train, ...
        'Learners', 'svm', ...
        'Coding', 'onevsall', ...  
        'Verbose', 1);       
    Y_pred = predict(svm_model, X_test);    
    accuracy = sum(Y_pred == Y_test) / length(Y_test);
    acc_wm_subtasks(fold,1) = accuracy;
end

%% Permutation Test
n_permutations = folds;
perm_accuracies = zeros(n_permutations, 1);
Y_true = Y;
for perm = 1:n_permutations
    Y_shuffled = Y_true(randperm(length(Y_true)));
    idx_train_perm = idx_train;
    idx_test_perm = idx_test;
    
    X_train_perm = X(idx_train_perm, :);
    Y_train_perm = Y_shuffled(idx_train_perm);
    X_test_perm = X(idx_test_perm, :);
    Y_test_perm = Y_shuffled(idx_test_perm);
    
    svm_model_perm = fitcecoc(X_train_perm, Y_train_perm, ...
        'Learners', 'svm', ...
        'Coding', 'onevsall', ...
        'Verbose', 0);  
    Y_pred_perm = predict(svm_model_perm, X_test_perm);
    perm_accuracies(perm) = sum(Y_pred_perm == Y_test_perm) / length(Y_test_perm);
    if mod(perm, 100) == 0
        fprintf('process: %d/%d\n', perm, n_permutations);
    end
end

% true_acuracy
true_accuracy = acc_wm_subtasks;
p_value = sum(perm_accuracies >= true_accuracy) / n_permutations;
figure;
histogram(perm_accuracies, 30, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
hold on;
xline(true_accuracy, 'r', 'LineWidth', 2);
xlabel('accuracy of calssification');
ylabel('count');
title(sprintf('Permutation Test Distribution (p = %.4f)', p_value));
legend({'Permutation Test Distribution', 'accuracy'});