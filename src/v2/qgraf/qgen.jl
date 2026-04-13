#  qgen — field assignment via dpntro lookup.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13797-14464
#
#  Phase 12a — dpntro builder (this file).
#  Phase 12b — vmap/lmap construction (extension to qg10).
#  Phase 12c — qgen recursive backtracker.

"""
    build_dpntro(vertex_rules) -> Dict{Int, Vector{Vector{Symbol}}}

Build the qgen lookup table.  Input is any iterable of `Vector{Symbol}`
representing expanded vertex rules (each Vector already sorted by
`_expand_vertex`'s `sort` step).  Output buckets the rules by arity and
sorts each bucket lexicographically for deterministic iteration.

Source: qgraf-4.0.6.f08:13889
  "vfo(vv) = stib(stib(dpntro(0)+vdeg(vv))+pmap(vv,1))"

The qgraf two-level nesting (degree → first-particle → list) is collapsed
to a single Dict-of-lists; downstream qgen filters at lookup time.  A
two-level dict is a future optimisation when profiling shows it matters.
"""
function build_dpntro(vertex_rules)
    dp = Dict{Int, Vector{Vector{Symbol}}}()
    for rule in vertex_rules
        d = length(rule)
        push!(get!(Vector{Vector{Symbol}}, dp, d), Vector{Symbol}(rule))
    end
    for d in keys(dp)
        sort!(dp[d])
    end
    dp
end

# Symmetric adjacency lookup from xg's stored upper triangle.
@inline _gam(state::TopoState, i::Integer, j::Integer) =
    i == j ? state.xg[i, i] : (i < j ? state.xg[i, j] : state.xg[j, i])

"""
    compute_qg10_labels(state) -> NamedTuple

Build the per-topology label arrays needed by qgen field assignment:
  vlis[i]   = vertex visited at position i (canonical traversal order)
  invlis[v] = inverse: position of vertex v in vlis
  rdeg[v]   = degree to already-visited vertices at the time v is picked
  sdeg[v]   = rdeg[v] + (self-loop edges at v)
  vmap[v,k] = the k-th neighbor of vertex v (sorted by invlis)
  lmap[v,k] = back-pointer: vmap[vmap[v,k], lmap[v,k]] == v

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12028-12102.

Externals (1..n_ext) are visited first in their natural order; internals
are visited greedily by maximum cumulative connection to the
already-visited set (qg10's `vaux` invariant).
"""
function compute_qg10_labels(state::TopoState)
    n     = Int(state.n)
    n_ext = Int(state.n_ext)
    rhop1 = Int(state.rhop1)

    vlis   = zeros(Int8, n)
    invlis = zeros(Int8, n)
    rdeg   = zeros(Int8, n)
    sdeg   = zeros(Int8, n)
    vmap   = zeros(Int8, n, n)
    lmap   = zeros(Int8, n, n)
    vaux   = zeros(Int8, n)

    # ── qg10:12028-12030 — vaux[i] = xn[i]
    @inbounds for i in 1:n
        vaux[i] = state.xn[i]
    end
    # qg10:12031-12034 — externals visited first, in their natural order.
    @inbounds for i in 1:n_ext
        vlis[i]   = Int8(i)
        invlis[i] = Int8(i)
    end
    # Internals start unvisited (invlis already zero).

    # ── qg10:12039-12069 — pick next internal by max vaux until all visited.
    aux = n_ext
    jj  = rhop1
    while aux < n
        # Skip already-visited.
        while jj <= n && invlis[jj] != Int8(0)
            jj += 1
        end
        # ii = argmax(vaux[k]) over unvisited k >= jj.
        ii = jj
        @inbounds for i in (jj + 1):n
            if invlis[i] == Int8(0) && vaux[i] > vaux[ii]
                ii = i
            end
        end
        vaux[ii] != Int8(0) || error("qg10_1 — no candidate vertex with positive vaux")
        aux += 1
        vlis[aux]   = Int8(ii)
        invlis[ii]  = Int8(aux)
        rdeg[ii]    = vaux[ii]
        sdeg[ii]    = rdeg[ii] + _gam(state, ii, ii)
        # Propagate: each remaining internal accumulates connections to ii.
        @inbounds for i in rhop1:n
            vaux[i] += _gam(state, i, ii)
        end
    end

    # ── qg10:12071-12102 — build vmap and lmap.
    fill!(vaux, Int8(0))
    @inbounds for i1 in 1:n
        ii = Int(vlis[i1])
        kk = 0
        j1 = 1
        while kk < Int(state.vdeg[ii])
            jj = Int(vlis[j1])
            sgn = Int8(1)
            for _ in 1:Int(_gam(state, ii, jj))
                kk += 1
                vmap[ii, kk] = Int8(jj)
                vaux[jj] += Int8(1)
                if ii != jj
                    lmap[ii, kk] = vaux[jj]
                else
                    lmap[ii, kk] = vaux[jj] + sgn
                    sgn = -sgn
                end
            end
            j1 += 1
        end
    end

    return (; vlis, invlis, vmap, lmap, rdeg, sdeg)
end
