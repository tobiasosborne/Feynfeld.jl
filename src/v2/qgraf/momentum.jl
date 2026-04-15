#  Momentum routing — spanning tree + edge classification.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13315-13558
#  Cross-ref: ALGORITHM.md §5.
#
#  Phase 16 minimal: spanning-tree builder + chord (loop-momentum-carrying)
#  identification.  Full momentum assignment via leaf-peeling deferred to
#  Layer-4 amplitude-evaluation work.

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

# ── Phase 18a-1: leaf-peel momentum routing ─────────────────────────────
#
# Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13400-13559
# Cross-ref: refs/qgraf/ALGORITHM.md §5.2

"""
    InternalEdge

One internal edge of a topology, with its routed momentum and edge type.

Fields:
- `v_lo, v_hi`: vertex endpoints (v_lo ≤ v_hi). v_lo == v_hi for self-loops.
- `parallel_idx`: 1..gam(v_lo, v_hi) for parallel edges.
- `momentum`: per-edge momentum (`nothing` if zero, single `Momentum` if
  one term with unit coeff, else `MomentumSum`).
- `edge_type`: `:rb` (regular bridge — externals only), `:sb` (special
  bridge — zero momentum), `:rnb` (regular non-bridge — chord with loop
  momentum), `:snb` (special non-bridge — self-loop chord).
"""
struct InternalEdge
    v_lo::Int
    v_hi::Int
    parallel_idx::Int
    momentum::Union{Nothing, Momentum, MomentumSum}
    edge_type::Symbol
end

"""
    EdgeMomenta

Per-topology momentum routing result.

Fields:
- `n_ext`: number of external legs (= length(ext_moms)).
- `ext_moms`: external momenta (mirrored input, qgraf "all incoming" convention).
- `internal`: list of `InternalEdge` for each internal edge in qgraf
  iteration order (vertex-pair (v_lo, v_hi) lexicographic, then by parallel index).
"""
struct EdgeMomenta
    n_ext::Int
    ext_moms::Vector{Momentum}
    internal::Vector{InternalEdge}
end

"""
    route_momenta(state, labels, ext_moms) -> EdgeMomenta

Per-edge momentum assignment via qgraf's leaf-peel algorithm.

# Convention
qgraf treats every external momentum as flowing INTO the diagram. The
caller is responsible for negating outgoing-physical momenta if a
physical-flow sign is desired.

# Inputs
- `state::TopoState` with adjacency `xg` and external counts `xn` filled.
- `labels` from `compute_qg10_labels(state)` (provides `vmap` for ext lookup).
- `ext_moms::Vector{Momentum}` of length `state.n_ext`.

# Algorithm (qgraf f08:13400-13559)
1. Init: each external edge i gets flow column i; each chord gets a fresh
   loop-momentum column; tree edges start at zero.
2. Leaf-peel: find leaf internal vertex (degree 1 in remaining tree),
   compute its tree edge as (sum of external contributions at the vertex)
   ± (sum of chord contributions, sign by endpoint match), remove the
   tree edge from the tree, restart.
3. Sign normalisation: if more than half of external coefficients on an
   edge are nonzero, flip per the f08:13511-13521 rule.
4. Edge-type tag: classify as :rb / :sb / :rnb / :snb per f08:13533-13558.

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13400-13559
Cross-ref: refs/qgraf/ALGORITHM.md §5.2
"""
function route_momenta(state::TopoState, labels,
                       ext_moms::Vector{Momentum};
                       loop_moms::Union{Nothing, Vector{Momentum}}=nothing)
    n_ext = Int(state.n_ext)
    length(ext_moms) == n_ext ||
        error("route_momenta: ext_moms length $(length(ext_moms)) ≠ state.n_ext $n_ext")
    n     = Int(state.n)
    rhop1 = Int(state.rhop1)

    # Enumerate internal edges in upper-triangular order. Self-loops contribute
    # `xg[i,i] ÷ 2` parallel copies (qgraf stores 2× the half-edge count).
    int_edges = Tuple{Int,Int,Int}[]
    int_first = Dict{Tuple{Int,Int}, Int}()
    for i in rhop1:n
        for j in i:n
            mult = i == j ? Int(state.xg[i,i]) ÷ 2 : Int(state.xg[i,j])
            mult == 0 && continue
            int_first[(i,j)] = length(int_edges) + 1
            for p in 1:mult
                push!(int_edges, (i, j, p))
            end
        end
    end
    nli_int = length(int_edges)
    nli_int == 0 && return EdgeMomenta(n_ext, ext_moms, InternalEdge[])

    # Phase 16 spanning tree picks the first parallel copy of each
    # tree-marked vertex pair; remaining parallels and self-loops are chords.
    tree_bool = build_spanning_tree(state)
    intree    = falses(nli_int)
    bb        = zeros(Int, n)
    for (k, (i, j, p)) in enumerate(int_edges)
        in_tree_k = (i != j) && tree_bool[i, j] && (p == 1)
        intree[k] = in_tree_k
        if in_tree_k
            bb[i] += 1
            bb[j] += 1
        end
    end

    # flow[k, c]: integer coefficient of momentum c on edge k.
    # cols 1..n_ext: external p_c; cols n_ext+1..end: loop momenta k_l.
    n_chord = nli_int - count(intree)
    flow    = zeros(Int, nli_int, n_ext + n_chord)
    cc = n_ext
    for k in 1:nli_int
        if !intree[k]
            cc += 1
            flow[k, cc] = 1
        end
    end

    # Resolve loop momentum names (default :k1, :k2, ...).
    loop_names = loop_moms === nothing ?
        [Momentum(Symbol(:k, l)) for l in 1:n_chord] :
        loop_moms
    length(loop_names) == n_chord ||
        error("route_momenta: loop_moms length $(length(loop_names)) ≠ n_chord $n_chord")

    _tail(k) = int_edges[k][1]
    _head(k) = int_edges[k][2]

    # Leaf-peel main loop (f08:13420-13491, labels 73, 71, 75).
    ii = rhop1
    @label loop_73
    if bb[ii] == 1
        j1 = 0
        # f08:13430-13437 — scan i1 < ii (lower-numbered internal neighbours).
        for i1 in rhop1:(ii-1)
            Int(state.xg[i1, ii]) > 0 || continue
            first = int_first[(i1, ii)]
            for k in first:(first + Int(state.xg[i1, ii]) - 1)
                if intree[k]; j1 = k; @goto found; end
            end
        end
        # f08:13438-13445 — scan i1 > ii.
        for i1 in (ii + 1):n
            Int(state.xg[ii, i1]) > 0 || continue
            first = int_first[(ii, i1)]
            for k in first:(first + Int(state.xg[ii, i1]) - 1)
                if intree[k]; j1 = k; @goto found; end
            end
        end
        error("route_momenta: leaf $ii has no tree edge (qg21_16)")
        @label found

        # f08:13451-13462 — external-momentum injection.
        if n_ext > 1 && Int(state.xn[ii]) > 0
            for i1 in 1:n_ext
                Int(labels.vmap[i1, 1]) == ii || continue
                if _tail(j1) == ii
                    flow[j1, i1] += 1
                else
                    flow[j1, i1] -= 1
                end
            end
        end

        # f08:13464-13486 — chord-momentum accumulation, sign by endpoint match.
        for i1 in rhop1:n
            i1 == ii && continue
            a, b = min(i1, ii), max(i1, ii)
            gij  = Int(state.xg[a, b])
            gij == 0 && continue
            first = get(int_first, (a, b), 0)
            first == 0 && continue
            for k in first:(first + gij - 1)
                intree[k] && continue
                if _head(j1) == _head(k) || _tail(j1) == _tail(k)
                    @inbounds for c in 1:size(flow, 2); flow[j1, c] -= flow[k, c]; end
                else
                    @inbounds for c in 1:size(flow, 2); flow[j1, c] += flow[k, c]; end
                end
            end
        end

        intree[j1]      = false
        bb[_tail(j1)]  -= 1
        bb[_head(j1)]  -= 1
        ii = rhop1
        @goto loop_73
    elseif ii < n
        ii += 1
        @goto loop_73
    end

    # Build InternalEdge results. Sign normalisation and full edge-type
    # classification arrive in subsequent TDD steps.
    internal = InternalEdge[]
    for k in 1:nli_int
        v_lo, v_hi, parallel_idx = int_edges[k]
        terms = Tuple{Rational{Int}, Momentum}[]
        for c in 1:n_ext
            flow[k, c] == 0 && continue
            push!(terms, (Rational{Int}(flow[k, c]), ext_moms[c]))
        end
        for l in 1:n_chord
            cc = n_ext + l
            flow[k, cc] == 0 && continue
            push!(terms, (Rational{Int}(flow[k, cc]), loop_names[l]))
        end
        mom = momentum_sum(terms)
        any_loop = any(flow[k, c] != 0 for c in (n_ext + 1):size(flow, 2))
        any_ext  = any(flow[k, c] != 0 for c in 1:n_ext)
        edge_type = any_loop ? (v_lo == v_hi ? :snb : :rnb) :
                    any_ext  ? :rb : :sb
        push!(internal, InternalEdge(v_lo, v_hi, parallel_idx, mom, edge_type))
    end
    EdgeMomenta(n_ext, ext_moms, internal)
end
