function state = dfs_exactL_sumabs(u, s, depth, sumAbsSoFar, L, adj, wabs, visited, path, edgeAbsPath, state)
    if depth == L
        % close u -> s
        nbrs = adj{u};
        wgts = wabs{u};
        for k = 1:numel(nbrs)
            if nbrs(k) == s
                closingAbs = wgts(k);
                score = sumAbsSoFar + closingAbs;

                edgeAbsFull = [edgeAbsPath(1:L-1), closingAbs];
                state = considerCycle_sumabs(path(1:L), edgeAbsFull, score, state);
                break
            end
        end
        return
    end

    nbrs = adj{u};
    wgts = wabs{u};

    for k = 1:numel(nbrs)
        v = nbrs(k);
        a = wgts(k);

        % enforce s is the smallest node in the cycle (avoid duplicates)
        if v < s, continue; end

        if ~visited(v)
            visited(v) = true;
            path(depth+1) = v;
            edgeAbsPath(depth) = a;
            state = dfs_exactL_sumabs(v, s, depth+1, sumAbsSoFar + a, L, adj, wabs, visited, path, edgeAbsPath, state);
            visited(v) = false;
        end
    end
end

