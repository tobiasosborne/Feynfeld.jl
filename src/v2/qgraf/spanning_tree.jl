#  Spanning tree + chord count (Phase 16 minimal).
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13315-13402
#  Cross-ref: ALGORITHM.md §5.

"""
    build_spanning_tree(state) -> Matrix{Bool}

Build a spanning tree over the topology's xg adjacency.  Returns a
Matrix{Bool} `in_tree[i, j]` (upper-triangular: i < j) indicating which
edges are in the tree.

Algorithm: BFS from vertex 1, picking edges greedily as encountered.
For each vertex, its first encounter via an edge contributes that edge
to the tree.  Subsequent edges to already-visited vertices are CHORDS.

Source: qgraf:13315-13402 (qg21 spanning-tree segment).  Note that qgraf
prefers high-multiplicity edges first; this minimal port uses plain BFS
order — equivalent for tree-vs-chord classification (the multiplicity
preference only affects which copy of a parallel edge is picked as tree).
"""
function build_spanning_tree(state::TopoState)
    n = Int(state.n)
    in_tree = falses(n, n)
    visited = falses(n)
    queue   = Int[]
    visited[1] = true
    push!(queue, 1)
    @inbounds while !isempty(queue)
        v = popfirst!(queue)
        for u in 1:n
            visited[u] && continue
            edge = u < v ? state.xg[u, v] : state.xg[v, u]
            edge > Int8(0) || continue
            in_tree[min(u, v), max(u, v)] = true
            visited[u] = true
            push!(queue, u)
        end
    end
    return in_tree
end

"""
    count_chords(state) -> Int

Number of chord edges — edges that close cycles (i.e., not in any
spanning tree).  Equals (P - V + 1) by Euler's formula = nloop.
"""
function count_chords(state::TopoState)
    n = Int(state.n)
    in_tree = build_spanning_tree(state)
    chord_count = 0
    @inbounds for i in 1:(n - 1)
        for j in (i + 1):n
            mult = Int(state.xg[i, j])
            mult == 0 && continue
            # `mult` parallel edges: 1 is in_tree (if visited), the rest are chords.
            if in_tree[i, j]
                chord_count += mult - 1
            else
                chord_count += mult
            end
        end
        # Self-loops are always chords (they form a 1-cycle).
        chord_count += Int(state.xg[i, i])  # (xg[i,i] = 2× self-loop count in qgraf)
    end
    chord_count += Int(state.xg[n, n])
    chord_count
end
