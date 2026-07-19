function state = considerCycle_sumabs(nodes, edgeAbsFull, score, state)
    % normalize by rotating so smallest node first
    [nodesN, edgeAbsN] = normalizeCycle_withEdges(nodes, edgeAbsFull);

    % de-dup by exact node sequence
    for t = 1:numel(state.score)
        if isequal(state.nodes{t}, nodesN)
            if score > state.score(t)
                state.nodes{t} = nodesN;
                state.score(t) = score;
                state.edgeAbs{t} = edgeAbsN;
            end
            return
        end
    end

    if numel(state.score) < state.K
        state.nodes{end+1} = nodesN;
        state.score(end+1) = score;
        state.edgeAbs{end+1} = edgeAbsN;
    else
        [minSc, minIdx] = min(state.score);
        if score > minSc
            state.nodes{minIdx} = nodesN;
            state.score(minIdx) = score;
            state.edgeAbs{minIdx} = edgeAbsN;
        end
    end
end
