function top = top_cycles_exactL_sumabs(W, L, K, tol)
% top_cycles_exactL_sumabs  Top-K strongest simple directed cycles of EXACT length L.
% Ranking score: sum of absolute edge weights along the cycle.
%
% S(c) = sum_{t=1..L} | w_{v_t, v_{t+1}} |, v_{L+1}=v_1
%
% Output 'top' fields:
%   .nodes  1×L nodes (rotated so smallest node first)
%   .len
%   .score  sum abs weights along the cycle
%   .edgeAbs  1×L abs weights along edges in order (including closing edge)

    if nargin < 4 || isempty(tol), tol = 0; end
    W = sparse(W);
    n = size(W,1);
    assert(size(W,2)==n, 'W must be square');
    assert(L>=2 && L==round(L), 'L must be an integer >=2');
    assert(K>=1 && K==round(K), 'K must be a positive integer');

    % Build adjacency list
    [ii,jj,vv] = find(W);
    keep = abs(vv) > tol;
    ii = ii(keep); jj = jj(keep); vv = vv(keep);

    adj = cell(n,1);
    wabs = cell(n,1);
    for e = 1:numel(ii)
        adj{ii(e)}(end+1) = jj(e);
        wabs{ii(e)}(end+1) = abs(vv(e));
    end

    % State for top-K
    state.K = K;
    state.nodes = cell(0,1);
    state.score = zeros(0,1);
    state.edgeAbs = cell(0,1);

    visited = false(n,1);
    path = zeros(1, L);
    edgeAbsPath = zeros(1, L-1); % for edges along path (excluding closing edge)

    for s = 1:n
        visited(:) = false;
        visited(s) = true;
        path(1) = s;
        state = dfs_exactL_sumabs(s, s, 1, 0, L, adj, wabs, visited, path, edgeAbsPath, state);
    end

    % Pack + sort
    m = numel(state.score);
    [~,ord] = sort(state.score, 'descend');

    top = struct('nodes', cell(m,1), 'len', cell(m,1), 'score', cell(m,1), 'edgeAbs', cell(m,1));
    for i = 1:m
        t = ord(i);
        top(i).nodes = state.nodes{t};
        top(i).len   = L;
        top(i).score = state.score(t);
        top(i).edgeAbs = state.edgeAbs{t};
    end
end

