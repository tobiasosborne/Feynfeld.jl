# Topological filters for Feynman diagram generation.
#
# Connectedness, 1PI (bridge-free), canonical form checking.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 3.6-3.7

# topology_types.jl defines FeynmanTopology, n_vertices

"""
    _is_connected_topo(topo) → Bool

Check if the topology represents a connected graph.
"""
function _is_connected_topo(topo::FeynmanTopology)
    n = n_vertices(topo)
    n <= 1 && return true
    # BFS from vertex 1
    visited = falses(n)
    queue = Int[1]
    visited[1] = true
    while !isempty(queue)
        v = popfirst!(queue)
        for u in 1:n
            if !visited[u] && (topo.adj[v, u] > 0 || topo.adj[u, v] > 0)
                visited[u] = true
                push!(queue, u)
            end
        end
    end
    all(visited)
end

"""
    _is_1pi(topo) → Bool

Check if the topology is 1-particle irreducible (no internal bridges).
A bridge is an internal edge whose removal disconnects the graph.
External edges (connecting external vertices to internal ones) are
always bridges and are NOT counted — 1PI refers to the internal structure.
"""
function _is_1pi(topo::FeynmanTopology)
    n = n_vertices(topo)
    n_ext = topo.n_ext
    # Only check INTERNAL edges (both endpoints are internal vertices)
    for i in (n_ext + 1):n, j in i:n
        topo.adj[i, j] > 0 || continue
        if topo.adj[i, j] == 1
            topo.adj[i, j] -= 1
            topo.adj[j, i] -= 1
            connected = _is_connected_topo(topo)
            topo.adj[i, j] += 1
            topo.adj[j, i] += 1
            connected || return false
        end
        # Multi-edges (multiplicity > 1) are never bridges
    end
    true
end

"""
    _is_canonical(adj, vdeg, n_ext) → Bool

Check if the adjacency matrix is in canonical form: the lexicographically
smallest representative of its isomorphism class.

For vertices in the same equivalence class, swap them and check if the
result is lexicographically smaller. If so, current form is NOT canonical.
"""
function _is_canonical(adj, vdeg, n_ext)
    n = length(vdeg)
    # For each pair of equivalent vertices, try swapping
    for i in 1:n, j in (i+1):n
        _same_class(vdeg, n_ext, i, j) || continue
        # Check if swapping i and j produces a lexicographically smaller matrix
        cmp = _compare_swapped(adj, n, i, j)
        cmp < 0 && return false  # swapped version is smaller → not canonical
    end
    true
end

"""
    _compare_swapped(adj, n, i, j) → Int

Compare adjacency matrix with its version where vertices i and j are swapped.
Returns -1 if swapped < original, 0 if equal, +1 if swapped > original.
"""
function _compare_swapped(adj, n, i, j)
    for r in 1:n, c in r:n
        # Map (r,c) through the swap permutation
        r2 = r == i ? j : r == j ? i : r
        c2 = c == i ? j : c == j ? i : c
        r2, c2 = minmax(r2, c2)  # keep upper triangle
        orig = adj[r, c]
        swapped = adj[r2, c2]
        orig < swapped && return 1   # original is smaller → keep
        orig > swapped && return -1  # swapped is smaller → not canonical
    end
    0  # identical (automorphism)
end

# Two vertices are in the same equivalence class if they have the same degree
# and are both internal (or both external).
function _same_class(vdeg, n_ext, i, j)
    both_ext = (i <= n_ext && j <= n_ext)
    both_int = (i > n_ext && j > n_ext)
    (both_ext || both_int) && vdeg[i] == vdeg[j]
end
