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

# ── Multiset helpers (simple Symbol → Int counter approach) ────────────

function _is_sub_multiset(sub::AbstractVector{Symbol}, sup::AbstractVector{Symbol})
    counts = Dict{Symbol, Int}()
    for f in sup
        counts[f] = get(counts, f, 0) + 1
    end
    for f in sub
        c = get(counts, f, 0)
        c == 0 && return false
        counts[f] = c - 1
    end
    true
end

function _multiset_diff(a::AbstractVector{Symbol}, b::AbstractVector{Symbol})
    counts = Dict{Symbol, Int}()
    for f in a
        counts[f] = get(counts, f, 0) + 1
    end
    for f in b
        counts[f] = get(counts, f, 0) - 1
    end
    result = Symbol[]
    for (f, c) in counts
        for _ in 1:c
            push!(result, f)
        end
    end
    sort!(result)
end

"""
    qgen_count_assignments(state, labels, ext_assignment, dpntro, conjugate) -> Int

Count valid field assignments to internal slots, given an external-slot
field assignment and the dpntro lookup table.  Multiset matching against
vertex rules; recursive backtracker over `labels.vlis` order.

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13880-13987.

Phase 12c minimal port:
  • Multiset rule matching (qgraf does positional matching with sorted
    rules; multiset is equivalent in algorithm, more idiomatic in Julia).
  • DOES NOT yet apply the symmetry-factor 1/S weighting that gives the
    true distinct-diagram count — that's Phase 15.  Returned count is the
    SUM over all (perm, assignment) tuples and over-counts by the
    topology automorphism size.
  • Self-loop slot validity (link(field) pairing) deferred until needed.
"""
function qgen_count_assignments(state::TopoState, labels,
                                 ext_assignment::AbstractVector{Symbol},
                                 dpntro::Dict{Int, Vector{Vector{Symbol}}},
                                 conjugate::Dict{Symbol, Symbol})
    n     = Int(state.n)
    n_ext = Int(state.n_ext)
    pmap  = fill(:_, n, MAX_V)

    # qgen:13880-13884 — externals: pmap[i,1] = ext field, propagate conjugate
    @inbounds for i in 1:n_ext
        f = ext_assignment[i]
        pmap[i, 1] = f
        nb   = Int(labels.vmap[i, 1])
        slot = Int(labels.lmap[i, 1])
        pmap[nb, slot] = conjugate[f]
    end

    return _qgen_recurse(state, labels, pmap, dpntro, conjugate, n_ext + 1)
end

function _qgen_recurse(state::TopoState, labels, pmap::Matrix{Symbol},
                       dpntro::Dict{Int, Vector{Vector{Symbol}}},
                       conjugate::Dict{Symbol, Symbol}, vind::Int)
    n = Int(state.n)
    if vind > n
        return 1   # all internal vertices satisfied
    end

    vv      = Int(labels.vlis[vind])
    deg     = Int(state.vdeg[vv])
    rdeg_vv = Int(labels.rdeg[vv])

    # Already-assigned fields at vv (slots 1..rdeg_vv from earlier vertices).
    assigned = Symbol[pmap[vv, k] for k in 1:rdeg_vv]
    sort!(assigned)

    rules = get(dpntro, deg, Vector{Vector{Symbol}}())
    count = 0

    # Save the slots we will mutate so we can restore on backtrack.
    saved = Vector{Tuple{Int, Int, Symbol}}()

    for rule in rules
        _is_sub_multiset(assigned, rule) || continue
        remaining = _multiset_diff(rule, assigned)

        # Assign remaining[1..end] to slots rdeg+1..deg (in given order),
        # propagate conjugate to neighbour at vmap[vv, slot] / lmap[vv, slot].
        empty!(saved)
        ok = true
        @inbounds for k in 1:length(remaining)
            slot = rdeg_vv + k
            pmap[vv, slot]  = remaining[k]
            nb   = Int(labels.vmap[vv, slot])
            nb_s = Int(labels.lmap[vv, slot])
            push!(saved, (nb, nb_s, pmap[nb, nb_s]))
            pmap[nb, nb_s] = conjugate[remaining[k]]
        end

        ok && (count += _qgen_recurse(state, labels, pmap, dpntro, conjugate, vind + 1))

        # Backtrack: restore neighbour slots; pmap[vv, rdeg+1..] are owned by us
        # and will be overwritten on next iteration (or left dangling — fine
        # since vind never revisits this vertex within this recursion frame).
        @inbounds for (nb, nb_s, prev) in saved
            pmap[nb, nb_s] = prev
        end
    end

    return count
end
