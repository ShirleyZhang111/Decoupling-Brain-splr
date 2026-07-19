function [nodesN, edgeAbsN] = normalizeCycle_withEdges(nodes, edgeAbsFull)
    % nodes: 1×L
    % edgeAbsFull: 1×L where edgeAbsFull(t)=|w_{nodes(t), nodes(t+1)}|, and nodes(L+1)=nodes(1)
    [~,p] = min(nodes);
    nodesN = [nodes(p:end), nodes(1:p-1)];
    edgeAbsN = [edgeAbsFull(p:end), edgeAbsFull(1:p-1)];
end