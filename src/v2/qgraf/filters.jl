#  qgraf filter predicates.
#
#  Source map (qgraf-4.0.6.f08):
#    qumpi(1,...)  3690  → onepi / nobridge
#    qumpi(2,...)         → nosbridge
#    qumpi(3,...)         → notadpole
#    qumpi(4,...)         → onshell
#    qumvi(1,...)  3777  → nosnail
#    qumvi(2,...)         → onevi
#    qumvi(3,...)         → onshellx
#    qgsig         13669 → nosigma
#    qcyc          18830 → cycli
#    inline in qg21       → noselfloop, nodiloop, noparallel
#
#  Phase 14a: the three inline xg-pattern filters (smallest, all in this file).

"""
    has_no_selfloop(state) -> Bool

True iff no internal vertex carries a self-loop (xg[i,i] == 0 for all
internal i).  Source: qg21:13960-13967 (dflag(dflaty%nosl) check).
"""
function has_no_selfloop(state::TopoState)
    @inbounds for i in Int(state.rhop1):Int(state.n)
        state.xg[i, i] > Int8(0) && return false
    end
    return true
end

"""
    has_no_diloop(state) -> Bool

True iff no two distinct vertices share more than one edge between them
(xg[i,j] ≤ 1 for all i < j).  Source: qg21:13969-13978 (dflag(dflaty%nodl)).
"""
function has_no_diloop(state::TopoState)
    n = Int(state.n)
    @inbounds for i in 1:(n - 1)
        for j in (i + 1):n
            state.xg[i, j] > Int8(1) && return false
        end
    end
    return true
end

"""
    has_no_parallel(state) -> Bool

True iff `has_no_diloop` AND no internal vertex has xg[i,i] > 3 (more
than one self-loop pair).  Source: qg21:13065-13076 (dflag(dflaty%nopa)).
"""
function has_no_parallel(state::TopoState)
    has_no_diloop(state) || return false
    @inbounds for i in Int(state.rhop1):Int(state.n)
        state.xg[i, i] > Int8(3) && return false
    end
    return true
end

# ── qumpi family — bridge filters (port of qg21:3690-3776) ─────────────

# Returns the set of vertices reachable from `start` in `state.xg`,
# treating xg as a symmetric adjacency, with edge (skip_i, skip_j)
# temporarily removed (set skip_i=0 to remove no edge).
function _bfs_reachable(state::TopoState, start::Int, skip_i::Int, skip_j::Int)
    n = Int(state.n)
    visited = falses(n)
    queue = Vector{Int}(undef, n)
    qhead = 1; qtail = 1
    visited[start] = true
    queue[qtail] = start; qtail += 1
    @inbounds while qhead < qtail
        v = queue[qhead]; qhead += 1
        for u in 1:n
            visited[u] && continue
            (v == skip_i && u == skip_j) && continue
            (v == skip_j && u == skip_i) && continue
            edge = u < v ? state.xg[u, v] : state.xg[v, u]
            edge > Int8(0) || continue
            visited[u] = true
            queue[qtail] = u; qtail += 1
        end
    end
    return visited
end

"""
    is_one_pi(state) -> Bool

True iff the topology is one-particle-irreducible: removing any single
internal edge does not disconnect the diagram.

Source: qgraf-4.0.6.f08:3690-3776 (qumpi(1, ...)).
"""
function is_one_pi(state::TopoState)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)
    @inbounds for i in rhop1:(n - 1)
        for j in (i + 1):n
            state.xg[i, j] == Int8(1) || continue
            visited = _bfs_reachable(state, 1, i, j)
            all(visited[1:n]) || return false
        end
    end
    return true
end

"""
    has_no_sbridge(state) -> Bool

True iff no internal-edge bridge has one component completely lacking
externals (= "self-bridge" in qgraf parlance).
Source: qg21:3690-3776 (qumpi(2, ...)).
"""
function has_no_sbridge(state::TopoState)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)
    n_ext = Int(state.n_ext)
    @inbounds for i in rhop1:(n - 1)
        for j in (i + 1):n
            state.xg[i, j] == Int8(1) || continue
            visited = _bfs_reachable(state, 1, i, j)
            all(visited[1:n]) && continue   # not a bridge
            ii = sum(visited[1:n_ext])      # externals in component containing 1
            (ii == 0 || ii == n_ext) && return false
        end
    end
    return true
end

"""
    has_no_tadpole(state) -> Bool

Same trigger as has_no_sbridge — qgraf uses xx=3 to RECORD the tadpole
without rejecting; the rejection variant matches xx=2.  Provided as a
distinct name for callers that want the notadpole semantic.
"""
has_no_tadpole(state::TopoState) = has_no_sbridge(state)

"""
    has_no_onshell(state) -> Bool

True iff no internal-edge bridge has one component with exactly 1 external.
Source: qg21:3690-3776 (qumpi(4, ...)).
"""
function has_no_onshell(state::TopoState)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)
    n_ext = Int(state.n_ext)
    @inbounds for i in rhop1:(n - 1)
        for j in (i + 1):n
            state.xg[i, j] == Int8(1) || continue
            visited = _bfs_reachable(state, 1, i, j)
            all(visited[1:n]) && continue
            ii = sum(visited[1:n_ext])
            (ii == 1 || ii == n_ext - 1) && return false
        end
    end
    return true
end

# ── qumvi family — vertex-cut filters (port of qg21:3777-3881) ─────────

"""
    has_no_snail(state) -> Bool

True iff no internal vertex carries a self-loop, with two qgraf-defined
exceptions where snails are LEGITIMATE (not filtered):
  • tree-level diagrams (nloop == 0): trivially no snails possible
  • 1-external 1-loop diagrams (rho(-1) == 1 AND nloop == 1):
    the snail IS the diagram (a legitimate self-energy insertion)

Source: qg21:3787-3804 (qumvi(1, ...)).
"""
function has_no_snail(state::TopoState)
    nloop = Int(state.nloop)
    nloop == 0 && return true
    n_ext = Int(state.n_ext)
    # qg21:3792-3797 — exception when rho(-1)==1 AND nloop==1.
    n_ext == 1 && nloop == 1 && return true
    @inbounds for i in Int(state.rhop1):Int(state.n)
        state.xg[i, i] != Int8(0) && return false
    end
    return true
end

# BFS that EXCLUDES a single vertex (used for vertex-cut detection).
function _bfs_skip_vertex(state::TopoState, start::Int, skip_v::Int)
    n = Int(state.n)
    visited = falses(n)
    queue   = Vector{Int}(undef, n)
    qhead   = 1; qtail = 1
    visited[start] = true
    visited[skip_v] = true     # treat as already-visited so it's never enqueued
    queue[qtail] = start; qtail += 1
    @inbounds while qhead < qtail
        v = queue[qhead]; qhead += 1
        for u in 1:n
            visited[u] && continue
            edge = u < v ? state.xg[u, v] : state.xg[v, u]
            edge > Int8(0) || continue
            visited[u] = true
            queue[qtail] = u; qtail += 1
        end
    end
    visited[skip_v] = false    # restore actual visited state for the skipped vertex
    return visited
end

"""
    is_one_vi(state) -> Bool

True iff the topology is one-vertex-irreducible: removing any single
INTERNAL vertex does not disconnect the diagram.

Source: qg21:3805-3881 (qumvi(2, ...)).
"""
function is_one_vi(state::TopoState)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)
    # qg21:3806-3808 — vacuous if there's at most one internal vertex.
    n - 1 - rhop1 <= 0 && return true
    @inbounds for vcut in rhop1:n
        # BFS from some other internal vertex (rhop1 if vcut≠rhop1 else n).
        start = (vcut != rhop1) ? rhop1 : n
        visited = _bfs_skip_vertex(state, start, vcut)
        # Check all OTHER internals reached (vcut itself is excluded).
        for i in rhop1:n
            i == vcut && continue
            visited[i] || return false   # vcut is a cut-vertex
        end
    end
    return true
end
