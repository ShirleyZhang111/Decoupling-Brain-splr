function [ICS,OCS,In_degree,Out_degree] = calc_node_wise_metrics(H)
% node-wise metrics calculation
% Input: 
%     H : Network connectivity matrix
% Output:
%     ICS: Incoming connection strength
%     OCS: Outgoing connection strength
%     In_degree: number of incoming connections
%     Out_degree: number of outgoing connections

[len,~] = size(H);
ICS = zeros(len,1);
OCS = zeros(len,1);
In_degree = zeros(len,1);
Out_degree = zeros(len,1);

for k = 1:len
    ICS(k,1) = sum(abs(H(k,:)));
    In_degree(k,1) = nnz(H(k,:));
    OCS(k,1) = sum(abs(H(:,k)));
    Out_degree(k,1) = nnz(H(:,k));
end
