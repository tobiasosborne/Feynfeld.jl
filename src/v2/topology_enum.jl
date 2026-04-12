# Topology enumeration via adjacency matrix construction.
#
# Builds all non-isomorphic connected multigraphs with prescribed
# vertex-degree sequence. Fills upper triangle of adjacency matrix
# entry by entry with backtracking.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 3

"""
    _enumerate_topologies(n_ext, partition; onepi=false) → Vector{FeynmanTopology}

Enumerate all non-isomorphic connected topologies with the given
external leg count and internal vertex-degree partition.
"""
function _enumerate_topologies(n_ext::Int, dp::DegreePartition; onepi::Bool=false)
    vdeg = Int[]
    for _ in 1:n_ext
        push!(vdeg, 1)
    end
    for (d, c) in sort(collect(dp.counts))
        for _ in 1:c
            push!(vdeg, d)
        end
    end
    n = length(vdeg)
    n == 0 && return FeynmanTopology[]

    results = FeynmanTopology[]
    adj = zeros(Int8, n, n)
    # Build list of (row, col) pairs in upper triangle to fill
    entries = Tuple{Int,Int}[]
    for i in 1:n, j in i:n
        push!(entries, (i, j))
    end
    _fill_entry!(results, adj, vdeg, n_ext, entries, 1, onepi)
    results
end

"""
    _fill_entry!(results, adj, vdeg, n_ext, entries, idx, onepi)

Backtracking: fill entry `idx` of the upper-triangle entry list.
Each entry (i,j) gets a value 0, 1, 2, ... up to the min of
remaining degrees of vertices i and j.
"""
function _fill_entry!(results, adj, vdeg, n_ext, entries, idx, onepi)
    n = length(vdeg)

    if idx > length(entries)
        # All entries filled — validate
        _validate_topology!(results, adj, vdeg, n_ext, onepi)
        return
    end

    i, j = entries[idx]
    rem_i = _remaining_degree(adj, vdeg, i, n)
    rem_j = _remaining_degree(adj, vdeg, j, n)

    if i == j
        # Self-loop: uses 2 half-edges from vertex i
        max_val = rem_i ÷ 2
    else
        max_val = min(rem_i, rem_j)
    end

    max_val = max(max_val, 0)

    for k in max_val:-1:0
        adj[i, j] = Int8(k)
        if i != j
            adj[j, i] = Int8(k)
        end
        _fill_entry!(results, adj, vdeg, n_ext, entries, idx + 1, onepi)
    end
    adj[i, j] = Int8(0)
    if i != j
        adj[j, i] = Int8(0)
    end
end

function _remaining_degree(adj, vdeg, v, n)
    used = 0
    for u in 1:n
        if u == v
            used += 2 * adj[v, v]
        else
            used += adj[v, u]
        end
    end
    vdeg[v] - used
end

function _validate_topology!(results, adj, vdeg, n_ext, onepi)
    n = length(vdeg)
    for v in 1:n
        _remaining_degree(adj, vdeg, v, n) == 0 || return
    end
    topo = FeynmanTopology(n_ext, copy(adj), copy(vdeg))
    _is_connected_topo(topo) || return
    onepi && !_is_1pi(topo) && return
    # Dedup via canonical form: only keep if this is the canonical representative
    _is_canonical_topo(topo) || return
    push!(results, topo)
end

"""
    _is_canonical_topo(topo) → Bool

Check if topo is the canonical (lex-smallest) representative of its
isomorphism class. Delegates to `QgrafPort.is_canonical_feynman`, which
runs the full per-equivalence-class lex-next-permutation check (qg21
labels 77/93/102/202/204 in qgraf-4.0.6.f08:13156-13291) replacing the
old pairwise-swap approach. Beads: feynfeld-ney, feynfeld-xjc.
"""
function _is_canonical_topo(topo::FeynmanTopology)
    QgrafPort.is_canonical_feynman(topo.adj, topo.vdeg, topo.n_ext)
end
