# Full per-equivalence-class lex-next-permutation canonicalisation.
#
# This is the bug locus: the old src/v2/topology_enum.jl does pairwise-swap
# canonicalisation which misses multi-vertex automorphisms, producing 474
# duplicate topologies at φ³ 2-loop φφ→φφ (qgraf says 465).
#
# Algorithm (Knuth Algorithm L applied per equivalence class in product):
#   1. Partition internal vertices into classes by (vdeg, xn, xg_diag).
#   2. Enumerate permutations that are products of class-internal permutations.
#      Advance rightmost class first; when exhausted, reset and advance next.
#   3. For each permutation π, compare gam(π(i), π(j)) vs gam(i, j) in
#      row-major order:
#        > 0 somewhere  →  current labelling NOT canonical (reject topology)
#        < 0 somewhere  →  π gives smaller image, continue search
#        = 0 everywhere →  π is an automorphism (ngsym += 1)
#   4. Externals (1..n_ext) are treated as distinguishable (each external is
#      a specific slot in the process) and never permuted. This matches
#      qgraf's effective behaviour after its external-leg fixity check.
#
# Ref: qgraf-4.0.6.f08:12819-12840 (xset/uset construction, label 19)
# Ref: qgraf-4.0.6.f08:13156-13291 (permutation + gam compare, labels 77/93/102/202/204/114/63)
# Ref: Nogueira1993_JCompPhys105_279.pdf §3 "Symmetry"
# Ref: ALGORITHM.md §3.5

"""
    compute_equiv_classes!(state) → Int

Partition internal vertices (rhop1..n) into equivalence classes by
(vdeg, xn, xg_diag). Does NOT assume sorted input: groups non-contiguous
same-invariant vertices into a shared class. Populates `state.classes`
and `state.xset`.

Returns the number of classes.
"""
function compute_equiv_classes!(state::TopoState)::Int
    empty!(state.classes)
    n = Int(state.n)
    rhop1 = Int(state.rhop1)
    rhop1 > n && return 0

    # Group internal vertices by (vdeg, xn, xg_diag).
    groups = Dict{Tuple{Int8,Int8,Int8}, Vector{Int8}}()
    for i in rhop1:n
        key = (state.vdeg[i], state.xn[i], state.xg[i, i])
        push!(get!(() -> Int8[], groups, key), Int8(i))
    end

    # Sort classes by key for deterministic order.
    for (key, members) in sort!(collect(groups); by = first)
        sort!(members)
        push!(state.classes, EquivClass(members, key[1]))
    end

    # Populate xset: class index per vertex.
    for (idx, cls) in enumerate(state.classes)
        for v in cls.members
            state.xset[Int(v)] = Int8(idx)
        end
    end
    length(state.classes)
end

"""
    _lex_next!(perm, first, last) → Bool

In-place lexicographic next-permutation on `perm[first..last]`
(Knuth Algorithm L, TAoCP §7.2.1.2). Returns `false` if `perm[first..last]`
was already the reverse-sorted (maximum) permutation.
"""
@inline function _lex_next!(perm::Vector{Int8}, first::Int, last::Int)::Bool
    first >= last && return false
    i = last - 1
    while i >= first && perm[i] >= perm[i + 1]
        i -= 1
    end
    i < first && return false
    j = last
    while perm[j] <= perm[i]
        j -= 1
    end
    perm[i], perm[j] = perm[j], perm[i]
    lo, hi = i + 1, last
    while lo < hi
        perm[lo], perm[hi] = perm[hi], perm[lo]
        lo += 1
        hi -= 1
    end
    true
end

"Reset perm within an equivalence class to the identity (members in order)."
@inline function _reset_identity!(perm::Vector{Int8}, cls::EquivClass)
    for (i, v) in enumerate(cls.members)
        perm[Int(v)] = v
    end
end

"""
    _lex_next_class!(perm, cls) → Bool

In-place lex-next-permutation of `perm` restricted to the (possibly
non-contiguous) vertex indices in `cls.members`. Reads/writes only at
those indices. Returns `false` when the class's permutation is the
reverse-sorted maximum.
"""
function _lex_next_class!(perm::Vector{Int8}, cls::EquivClass)::Bool
    mem = cls.members
    k = length(mem)
    k <= 1 && return false
    # Gather current values at the class's positions.
    # Find the rightmost i such that perm[mem[i]] < perm[mem[i+1]].
    i = k - 1
    while i >= 1 && perm[Int(mem[i])] >= perm[Int(mem[i + 1])]
        i -= 1
    end
    i < 1 && return false
    # Find rightmost j > i with perm[mem[j]] > perm[mem[i]].
    j = k
    while perm[Int(mem[j])] <= perm[Int(mem[i])]
        j -= 1
    end
    # Swap perm[mem[i]] <-> perm[mem[j]].
    a = Int(mem[i]); b = Int(mem[j])
    perm[a], perm[b] = perm[b], perm[a]
    # Reverse perm at positions mem[i+1..k].
    lo, hi = i + 1, k
    while lo < hi
        pa = Int(mem[lo]); pb = Int(mem[hi])
        perm[pa], perm[pb] = perm[pb], perm[pa]
        lo += 1
        hi -= 1
    end
    true
end

"""
    next_class_perm!(perm, classes) → Bool

Product-space lex-next-permutation across all equivalence classes. Advances
the rightmost class whose in-class permutation can still advance, resets
all later classes to identity. Returns `false` when the entire product
space is exhausted.
"""
function next_class_perm!(perm::Vector{Int8}, classes::Vector{EquivClass})::Bool
    for k in length(classes):-1:1
        cls = classes[k]
        if _lex_next_class!(perm, cls)
            for j in (k + 1):length(classes)
                _reset_identity!(perm, classes[j])
            end
            return true
        end
        _reset_identity!(perm, cls)
    end
    false
end

"""
    _compare_permuted_adjacency(state, perm) → Int

Compare `gam(π(i), π(j))` vs `gam(i, j)` in row-major order for internal
vertex pairs (i ≥ rhop1, j ≥ i). Returns `+1` on first strict increase,
`-1` on first strict decrease, `0` if all equal.

qg21 stores only the upper triangle in `xg`; `gam` is the symmetric
adjacency, accessed via `xg[min, max]`.
"""
function _compare_permuted_adjacency(state::TopoState, perm::Vector{Int8})::Int
    n = Int(state.n)
    xg = state.xg
    # Iterate over ALL vertex pairs (externals included). Externals are fixed
    # by perm (perm[i]=i for i ≤ n_ext), so ext-ext entries compare trivially.
    # Ext-int entries are what distinguish labelings of the same physical
    # topology under internal-vertex permutation — they MUST be compared.
    for i in 1:n
        pi_ = Int(perm[i])
        for j in i:n
            pj_ = Int(perm[j])
            g_perm = pi_ <= pj_ ? xg[pi_, pj_] : xg[pj_, pi_]
            g_orig = xg[i, j]
            if g_perm > g_orig
                return +1
            elseif g_perm < g_orig
                return -1
            end
        end
    end
    0
end

"""
    is_canonical_full!(state) → Bool

Returns `true` iff the current labelling of `state.xg` is the lex-smallest
representative of its isomorphism class under permutations of equivalence
classes. Side effect: sets `state.ngsym` to the automorphism count.

This replaces the old `_is_canonical_topo` (pairwise swaps only, bug locus
for the 474/465 φ³ 2-loop overcount).
"""
function is_canonical_full!(state::TopoState)::Bool
    compute_equiv_classes!(state)
    perm = state.perm_buf
    for i in 1:Int(state.n)
        perm[i] = Int8(i)
    end
    state.ngsym = Int32(1)  # identity is always an automorphism
    while next_class_perm!(perm, state.classes)
        cmp = _compare_permuted_adjacency(state, perm)
        # Convention: canonical = lex-SMALLEST representative (matches
        # topology_enum.jl's descending fill-and-accept iteration).
        # cmp > 0 means permuted > orig — keep searching, orig still may win.
        # cmp < 0 means permuted < orig — a smaller relabelling exists → reject.
        # cmp = 0 means automorphism.
        if cmp < 0
            return false
        elseif cmp == 0
            state.ngsym += Int32(1)
        end
    end
    true
end

"""
    is_canonical_feynman(adj, vdeg, n_ext) → Bool

Bridge: canonicality check for the legacy `FeynmanTopology` shape produced
by `topology_enum.jl` (which fills adjacency entries without qg21's Step
A/B sorting). Because vertices within a degree are NOT pre-sorted by
(xn, xg_diag), we use DEGREE-ONLY equivalence classes and let the full
per-class lex-next-permutation check discover the canonical ordering.

This is a drop-in replacement for the old pairwise-swap `_is_canonical_topo`.
The full-permutation expansion is what fixes the 474→465 bug at φ³ 2-loop:
where pairwise swap missed multi-vertex automorphisms, full S_n enumeration
catches them.

One TopoState is allocated per call; hot-path uses go through
`is_canonical_full!(state)` directly with pre-computed classes.
"""
function is_canonical_feynman(adj::AbstractMatrix{<:Integer},
                               vdeg::AbstractVector{<:Integer},
                               n_ext::Integer)::Bool
    n = length(vdeg)
    n ≤ MAX_V || error("is_canonical_feynman: n=$n exceeds MAX_V=$MAX_V")
    n_int = n - n_ext
    n_int <= 0 && return true   # no internals to permute
    p = Partition(Int8(n_ext), Int8[max(n_int, 0)], Int8(3), Int8(0))
    s = TopoState(p)
    s.n = Int8(n)
    s.n_ext = Int8(n_ext)
    s.rhop1 = Int8(n_ext + 1)
    fill!(s.vdeg, Int8(0))
    for i in 1:n
        s.vdeg[i] = Int8(vdeg[i])
    end
    fill!(s.xg, Int8(0))
    for i in 1:n, j in i:n
        s.xg[i, j] = Int8(adj[i, j])
    end
    # Build DEGREE-ONLY classes. Internals are already degree-sorted by
    # topology_enum.jl construction, so contiguous by degree.
    empty!(s.classes)
    rhop1 = Int(s.rhop1)
    first_v = Int8(rhop1)
    cur_deg = s.vdeg[rhop1]
    for i in (rhop1 + 1):n
        if s.vdeg[i] != cur_deg
            push!(s.classes, EquivClass(first_v, Int8(i - 1), cur_deg))
            first_v = Int8(i)
            cur_deg = s.vdeg[i]
        end
    end
    push!(s.classes, EquivClass(first_v, Int8(n), cur_deg))

    # Run the full permutation check with these classes.
    perm = s.perm_buf
    for i in 1:n
        perm[i] = Int8(i)
    end
    s.ngsym = Int32(1)
    while next_class_perm!(perm, s.classes)
        cmp = _compare_permuted_adjacency(s, perm)
        if cmp < 0
            return false
        elseif cmp == 0
            s.ngsym += Int32(1)
        end
    end
    true
end
