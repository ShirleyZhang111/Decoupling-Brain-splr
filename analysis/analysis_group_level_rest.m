% Group-level analysis of resting-state coupling components (S+rank-1)
load('../results/group_rest_Coup_rank1.mat')
[n,rank,items] = size(U_rest);

%% Analysis of Generative low-rank component
% U: region * rank * subject
[coeff_u, score_u, latent_u, tsquared_u, explained_u] = pca(reshape(U_rest,379,[])');
figure;
plot(explained_u(1:10),'.-','Color',[0,0,0],'LineWidth',1.5,'MarkerSize',5);
xlim([0.5,10.5]);
ylim([-2,85]);
% group-averaged u_1
u1 = coeff_u(:,1); 
% group-averaged u_2
u2 = coeff_u(:,2);
% group-averaged u_3
u3 = coeff_u(:,3);

% V_hat_all: region * subject
% Group-averaged u is significantly correlated with 
% the group-averaged real part of embedded v_hat.
% [coeff_v_hat, score_v_hat, latent_v_hat, tsquared_v_hat, explained_v_hat] = pca(V_hat_all.');
% [r1,p1] = corr(coeff_u(:,1),coeff_v_hat(:,1));

%% Analysis of Generative sparse component (Figure.3)
load('../data/7Networks.mat')
H = zeros(n,n);
S = zeros(n,n);
for k = 1:items
    H = H + Sparse_Coup_rest(:,:,k) + U_rest(:,k)*V_rest(:,k)';
    S = S + Sparse_Coup_rest(:,:,k);
end
H_mean = H/k;
S_mean = S/k;
list = [Visual;Somatomotor;Ventral;Dorsal;Limbic;Frontal;DMN];

% Generative connectivity (S + u*v')
figure;
imagesc(H_mean(list,list));
axis square;
clim([-0.003,0.003])
colormap(slanCM('bwr'));
setPivot(0);

% Generative sparse component (S)
figure;
imagesc(S_mean(list,list));
axis square;
clim([-0.003,0.003])
colormap(slanCM('bwr'));
setPivot(0);

% Modular organization of the generative sparse component.
load('../data/HCPex_LabelID.mat');
modular = sortrows(cell2mat(LabelID(:,[2,5])),1);

Modu_pos = zeros(items,1);
Modu_neg = zeros(items,1);
[m,~,items] = size(Sparse_Coup);
for id = 1:items
    Coup_pos = Sparse_Coup(:,:,id)- diag(diag(Sparse_Coup(:,:,id)));
    Coup_pos(Coup_pos < 0) = 0;
    [Modu_pos(id,1),~] = calc_network_metrics(Coup_pos,modular);
    
    Coup_neg = Sparse_Coup(:,:,id)- diag(diag(Sparse_Coup(:,:,id)));
    Coup_neg(Coup_neg > 0) = 0;
    [Modu_neg(id,1),~] = calc_network_metrics(Coup_neg,modular);
end

all_data = [Modu_neg; Modu_pos];
group = [ones(length(Modu_neg),1); 2*ones(length(Modu_pos),1)];

figure;
hb = boxplot(all_data,group,...              
                    'Color','k',...                                   
                    'symbol',' ',...                                  
                    'Notch','on',...                                  
                    'OutlierSize',4);                               

colors = [151 208 197;151 208 197]/255;
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),colors(j,:),'FaceAlpha',.9);
end
box on
set(gca, 'Box', 'off');