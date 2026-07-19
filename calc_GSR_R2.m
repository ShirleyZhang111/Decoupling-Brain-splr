function  R2_GSR2 = calc_GSR_R2(inCfg)

data_path = inCfg.data_path;
file_name = inCfg.file_name;
if isfield(inCfg,'T');T = inCfg.T;else;T = [];end

% start processing
load([data_path file_name],'TC');% N*T
TC = TC(:,6:end-5);
if isempty(T);T = 5:size(TC,2)-4;end
data = TC(:,T)';
[T,N] = size(data);

% global signal: g
g = mean(data, 2); % [T x 1]
residuals = zeros(T, N); 

for i = 1:N
    yi = data(:, i);
    
    % linear regression ：yi = beta_i * g + epsilon
    beta_i = (g' * g) \ (g' * yi); 
    y_pred = beta_i * g;
    epsilon_i = yi - beta_i * g;
  
    beta(i,1) = beta_i;
    residuals(:, i) = epsilon_i;
    
    var_yi = var(yi, 1); % 总方差
    var_epsilon = var(epsilon_i, 1); % 残差方差
    R2_GSR2(i,1) = 1 - var_epsilon / var_yi;
end
