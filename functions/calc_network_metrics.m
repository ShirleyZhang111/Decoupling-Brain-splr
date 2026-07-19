function [Q,r] = calc_network_metrics(H,modular)
% Network metrics calculation
% Input: 
%     H : Network connectivity matrix
%     modular : Module distribution of the network
% Output:
%     Q : Modularity Index
%     r : Assortativity

[len,~] = size(H);
delta = zeros(len,len);
for i = 1:len
    for j = 1:len
        delta(i,j) = modular(i,2) == modular(j,2);
    end
end
delta = double(delta);

E_dir = sum(sum(abs(H))); 
A = abs(H);
k = sum(A,2);

% Modularity
q = zeros(len,len);
for i = 1:len
    for j = 1:len
        q(i,j) = (A(i,j)-k(i)*k(j)/E_dir)*delta(i,j);
    end
end
Q = sum(sum(q))/E_dir/2;

% Assortativity
r1 = 0; r2 = 0; r3 = 0; 
for i = 1:379
    for j = i:379
        r1 = r1 + k(i)*k(j)*A(i,j);
        r2 = r2 + 0.5*(k(i)+k(j))*A(i,j);
        r3 = r3 + 0.5*(k(i)^2+k(j)^2)*A(i,j);
    end
end
r = (r1/E_dir - (r2/E_dir)^2)/(r3/E_dir- (r2/E_dir)^2);
