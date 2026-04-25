#  Phase 17 dedup audition: three strategies, same emission stream.
#
#  For phi3 φφ → φφ tree, qg21+qg10+qgen produces 24 raw emissions but
#  qgraf reports 3 distinct diagrams (s, t, u channels).  The factor-8
#  collapse comes from the topology automorphism group acting on the
#  external-leg labelling.
#
#  qgraf's actual approach (per qgen:14036-14100) is approach (B):
#  rejection-based canonical-pmap check.  For Phase 17 audition we
#  implement all three:
#
#    (A) Burnside       Σ 1/|Stab(emission)| over emissions
#    (B) Canonical-pmap emit only when no auto gives a lex-smaller diagram
#    (C) Pre-filter     enumerate only one perm per orbit at qg10 stage
#
#  All three should agree on every diagram count.  (B) is the qgraf-
#  faithful path; (A) gives a Burnside cross-check; (C) is the fastest
#  in principle (no per-emission auto check).

using Combinatorics: permutations as _permutations

# ── Per-vertex signature: sorted (perm-applied-neighbour, field) pairs ──
function _vertex_sig(state::TopoState, labels, perm::AbstractVector{<:Integer},
                     pmap::AbstractMatrix{Symbol}, v_orig::Int)
    vd = Int(state.vdeg[v_orig])
    pairs = Vector{Tuple{Int8, Symbol}}(undef, vd)
    @inbounds for slot in 1:vd
        nb = Int(labels.vmap[v_orig, slot])
        pairs[slot] = (Int8(perm[nb]), pmap[v_orig, slot])
    end
    sort!(pairs)
    pairs
end

# Diagram signature: sigs[k] = signature at NEW vertex k = π(v_orig).
function _diagram_sig(state::TopoState, labels, perm::AbstractVector{<:Integer},
                      pmap::AbstractMatrix{Symbol})
    n = Int(state.n)
    inv_perm = zeros(Int8, n)
    @inbounds for v in 1:n
        inv_perm[Int(perm[v])] = Int8(v)
    end
    sigs = Vector{Vector{Tuple{Int8, Symbol}}}(undef, n)
    @inbounds for k in 1:n
        sigs[k] = _vertex_sig(state, labels, perm, pmap, Int(inv_perm[k]))
    end
    sigs
end

# Inverse of a vertex permutation — used by every action below to map
# new-vertex coordinates back to old-vertex coordinates.
function _inv_perm(perm::AbstractVector{<:Integer})
    n = length(perm)
    out = Vector{Int8}(undef, n)
    @inbounds for v in 1:n
        out[Int(perm[v])] = Int8(v)
    end
    out
end

# ── Right-action convention on (ps1, pmap) ────────────────────────────────
#
# A topology automorphism g ∈ S_n permutes vertices.  External slots
# 1..n_ext are themselves vertices, so g induces a permutation of slot
# indices.  The PHYSICAL leg identification (which index is in_e, etc.)
# does NOT change under g — autos describe combinatorial symmetry, not
# physical relabeling.
#
# Therefore the action of g on `ps1 : slot → physical-leg` is:
#     (g·ps1)[i] = ps1[g⁻¹[i]]               (right action)
# i.e. "the new slot i was old slot g⁻¹[i], whose physical leg was
# ps1[g⁻¹[i]]".  The action on `pmap` is consistent: `_diagram_sig`
# already views pmap at new vertex k via old vertex inv_perm[k].
#
# The earlier left-action `(g·ps1)[i] = g[ps1[i]]` is a different (and
# physically wrong) action — it relabels physical legs, not slots.  Under
# left action, separate qgen emissions in the same physical orbit (e.g.
# Bhabha t-channel members connected by topology auto (3,4)) are NOT
# identified, and `is_emission_canonical` over-accepts.
#
# Verified empirically: with right action, is_emission_canonical accepts
# exactly the qgraf-validated diagram count (2 for Bhabha, 1 for ee→μμ,
# matching count_diagrams_qg21).  See bead feynfeld-vjw9.

"""
    _ps1_preserved(perm, ps1) -> Bool

True iff `perm` stabilises `ps1` under the right action
`(perm·ps1)[i] = ps1[perm⁻¹[i]]`.  Body checks the equivalent condition
`ps1[perm[i]] == ps1[i]` for all `i` — this is `ps1 ∘ perm == ps1`
(right-composing both sides of `ps1 ∘ perm⁻¹ == ps1` by `perm`).  For a
bijective `ps1` this collapses to `perm = identity` on the external slots.
"""
function _ps1_preserved(perm::AbstractVector{<:Integer},
                         ps1::AbstractVector{<:Integer})
    @inbounds for i in eachindex(ps1)
        Int(ps1[Int(perm[i])]) == Int(ps1[i]) || return false
    end
    true
end

"""
    is_emission_canonical(state, labels, autos, ps1, pmap) -> Bool

True iff the original labelling (identity perm) gives the LEX-SMALLEST
`(ps1, diagram-signature)` pair under the right action of `autos` —
i.e. for every non-identity `g ∈ autos`, the pair
`(ps1∘g⁻¹, _diagram_sig(g, pmap))` is ≥ `(ps1, _diagram_sig(id, pmap))`.

This is qgraf's full dedup criterion (qgen:14036-14100), generalised to
account for the external-leg labelling that ps1 encodes.  Each orbit
contains exactly one canonical emission, so for any process the count
of accepted emissions equals `count_diagrams_qg21(...)`.
"""
function is_emission_canonical(state::TopoState, labels,
                                autos::Vector{Vector{Int8}},
                                ps1::AbstractVector{<:Integer},
                                pmap::AbstractMatrix{Symbol})
    n = Int(state.n)
    identity = Int8.(1:n)
    orig_pmap_sig = _diagram_sig(state, labels, identity, pmap)
    n_ext = length(ps1)
    orig = (Vector{Int}(ps1), orig_pmap_sig)
    @inbounds for i in 2:length(autos)
        auto = autos[i]
        ginv = _inv_perm(auto)
        new_ps1 = Int[Int(ps1[Int(ginv[j])]) for j in 1:n_ext]
        new_pmap_sig = _diagram_sig(state, labels, auto, pmap)
        new_pair = (new_ps1, new_pmap_sig)
        new_pair < orig && return false
    end
    true
end

"""
    same_emission_orbit(state, labels, autos, ps1_a, pmap_a, ps1_b, pmap_b) -> Bool

True iff there exists an automorphism `g ∈ autos` mapping emission A to
emission B under the right action: `ps1_a ∘ g⁻¹ == ps1_b` AND
`_diagram_sig(g, pmap_a) == _diagram_sig(id, pmap_b)`.  Uses the same
`_diagram_sig` view that `emission_stabilizer` and `is_emission_canonical`
use, so all three agree on the orbit equivalence relation.

Cost: O(|autos| · (n_ext + Σ vdeg)).  For tree 2→2 problems, |autos|
and emission count are both small.
"""
function same_emission_orbit(state::TopoState, labels,
                               autos::Vector{Vector{Int8}},
                               ps1_a::AbstractVector{<:Integer},
                               pmap_a::AbstractMatrix{Symbol},
                               ps1_b::AbstractVector{<:Integer},
                               pmap_b::AbstractMatrix{Symbol})
    length(ps1_a) == length(ps1_b) || return false
    n_ext = length(ps1_a)
    n     = Int(state.n)
    identity = Int8.(1:n)
    target_sig = _diagram_sig(state, labels, identity, pmap_b)
    @inbounds for g in autos
        ginv = _inv_perm(g)
        # ps1 right-action: ps1_a ∘ g⁻¹ == ps1_b ?
        ps1_ok = true
        for j in 1:n_ext
            if Int(ps1_a[Int(ginv[j])]) != Int(ps1_b[j]); ps1_ok = false; break; end
        end
        ps1_ok || continue
        # pmap action via _diagram_sig: g · pmap_a ≡ pmap_b ?
        if _diagram_sig(state, labels, g, pmap_a) == target_sig
            return true
        end
    end
    false
end

"""
    emission_stabilizer(state, labels, autos, ps1, pmap) -> Int

Number of autos that preserve BOTH the pmap signature AND the ps1
labelling — i.e., |Stab(emission)| in the orbit-stabilizer sense over the
(ps1, pmap) joint space.  Identity always counts.
"""
function emission_stabilizer(state::TopoState, labels,
                              autos::Vector{Vector{Int8}},
                              ps1::AbstractVector{<:Integer},
                              pmap::AbstractMatrix{Symbol})
    n = Int(state.n)
    identity = Int8.(1:n)
    orig_pmap_sig = _diagram_sig(state, labels, identity, pmap)
    n_inv = 1
    @inbounds for i in 2:length(autos)
        auto = autos[i]
        _ps1_preserved(auto, ps1) || continue
        new_pmap_sig = _diagram_sig(state, labels, auto, pmap)
        new_pmap_sig == orig_pmap_sig && (n_inv += 1)
    end
    n_inv
end

# Backwards-compat aliases (older name; same semantics as the *_canonical
# variants without ps1 — used only by the failing first-cut audition).
is_pmap_canonical(state, labels, autos, pmap) =
    is_emission_canonical(state, labels, autos, 1:Int(state.n_ext), pmap)
pmap_stabilizer(state, labels, autos, pmap) =
    emission_stabilizer(state, labels, autos, 1:Int(state.n_ext), pmap)

# ── Common emission stream: yields (state, labels, ps1_perm, pmap) ─────
#
# Uses the legacy partition iterator (Phase 10 not yet ported), then runs
# qg21+qg10+qgen, calling `callback(state, labels, ps1, pmap)` for each
# valid emission.  pmap is mutated in place by qgen — copy if needed.
#
# `ext_exp::Union{Nothing, Vector{Symbol}}` lets callers supply the
# per-slot field labels (qgraf "all incoming" convention) directly,
# bypassing `_expand_external_fields`. Required when the alternation
# heuristic would mis-label legs (e.g. Bhabha — see
# `Feynfeld._qgraf_ext_labels`). When `nothing`, falls back to the
# alternation expansion of `vcat(in_fields, out_fields)`.
function _foreach_emission(callback::F, model, in_fields::Vector{Symbol},
                            out_fields::Vector{Symbol}; loops::Int,
                            ext_exp::Union{Nothing, Vector{Symbol}}=nothing) where {F}
    n_ext   = length(in_fields) + length(out_fields)
    ext_raw = vcat(in_fields, out_fields)
    exp     = _expand_model_for_diagen(model)
    ext_exp_eff = ext_exp === nothing ? _expand_external_fields(ext_raw, exp) :
                                          ext_exp
    length(ext_exp_eff) == n_ext ||
        error("_foreach_emission: expanded labels length $(length(ext_exp_eff)) ≠ n_ext $n_ext")
    dpntro  = build_dpntro(exp.vertex_rules)
    rules   = feynman_rules(model)
    vd      = Set(length(k) for k in keys(rules.vertices))
    for dp in _degree_partitions(n_ext, loops, vd)
        degs = sort([d for (d, c) in dp.counts if c > 0])
        if isempty(degs); continue; end
        mrho = degs[1]
        cv   = Int8[get(dp.counts, d, 0) for d in mrho:degs[end]]
        qp   = Partition(Int8(n_ext), cv, Int8(mrho), Int8(loops))
        s    = TopoState(qp)
        qg21_enumerate!(s) do state
            labels = compute_qg10_labels(state)
            ps1 = collect(1:n_ext)
            while true
                ext_perm = Symbol[ext_exp_eff[ps1[i]] for i in 1:n_ext]
                qgen_enumerate_assignments(state, labels, ext_perm, dpntro,
                                            exp.conjugate) do st, pmap
                    callback(st, labels, ps1, pmap)
                end
                # next perm
                j = n_ext - 1
                while j >= 1 && ps1[j] >= ps1[j + 1]; j -= 1; end
                j == 0 && break
                k = n_ext
                while ps1[k] <= ps1[j]; k -= 1; end
                ps1[j], ps1[k] = ps1[k], ps1[j]
                reverse!(view(ps1, (j + 1):n_ext))
            end
        end
    end
end

"""
    count_dedup_burnside(model, in_fields, out_fields; loops) -> Int

Approach (A) — Burnside: Σ 1/|Stab(emission)| over all emissions.  The sum
is exact when emissions partition into orbits (which they do under topology
auto action), giving #orbits = distinct-diagram count.
"""
function count_dedup_burnside(model, in_fields::Vector{Symbol},
                                out_fields::Vector{Symbol}; loops::Int)
    # #orbits = (1/|G|) × Σ_e |Stab(e)|
    # where the sum is over all (ps1, pmap) emissions and G is the
    # topology automorphism group acting on (ps1, pmap) jointly.
    total = Rational{Int}(0)
    _foreach_emission(model, in_fields, out_fields; loops=loops) do state, labels, ps1, pmap
        autos = enumerate_topology_automorphisms(state)
        g_size = length(autos)
        n_stab = emission_stabilizer(state, labels, autos, ps1, pmap)
        total += n_stab // g_size
    end
    @assert denominator(total) == 1 "Burnside sum non-integer ($total)"
    Int(numerator(total))
end

"""
    count_dedup_canonical(model, in_fields, out_fields; loops) -> Int

Approach (B) — qgraf-faithful canonical-pmap rejection: count emissions
whose diagram signature is the lex-smallest in its automorphism orbit.
"""
function count_dedup_canonical(model, in_fields::Vector{Symbol},
                                 out_fields::Vector{Symbol}; loops::Int)
    n = 0
    _foreach_emission(model, in_fields, out_fields; loops=loops) do state, labels, ps1, pmap
        autos = enumerate_topology_automorphisms(state)
        if is_emission_canonical(state, labels, autos, ps1, pmap)
            n += 1
        end
    end
    n
end

"""
    count_diagrams_qg21(model, in_fields, out_fields; loops, onepi=false,
                        nosbridge=false, notadpole=false, onshell=false,
                        nosnail=false, onevi=false,
                        noselfloop=false, nodiloop=false, noparallel=false) -> Int

Primary entry point for the qg21 Strategy C diagram counter.  Wraps
qg21 → qg10 → qgen and applies (A) Burnside dedup over (ps1, pmap) joint
orbits under the full topology auto group Γ_F.

Filter kwargs mirror qgraf's dflag options; each rejects topologies
failing the named predicate.  Filters not yet ported: `nosigma`,
`cycli`, `onshellx`, `floop` (fermion-loop required), `bipart`.

Ref: qg21 filter map (ALGORITHM.md §3.7).
"""
function count_diagrams_qg21(model, in_fields::Vector{Symbol},
                              out_fields::Vector{Symbol};
                              loops::Int=0,
                              onepi::Bool=false,
                              nosbridge::Bool=false,
                              notadpole::Bool=false,
                              onshell::Bool=false,
                              nosnail::Bool=false,
                              onevi::Bool=false,
                              noselfloop::Bool=false,
                              nodiloop::Bool=false,
                              noparallel::Bool=false)
    n_ext   = length(in_fields) + length(out_fields)
    ext_raw = vcat(in_fields, out_fields)
    exp     = _expand_model_for_diagen(model)
    # OK to use the alternation expansion here: this counter only computes
    # diagram cardinality, which depends on the field-label MULTISET (not on
    # which slot hosts which conjugate-pair member).  See gpi5 / Step 2 for
    # the pipeline-side fix where slot-physical alignment matters.
    ext_exp = _expand_external_fields(ext_raw, exp)
    dpntro  = build_dpntro(exp.vertex_rules)
    rules   = feynman_rules(model)
    vd      = Set(length(k) for k in keys(rules.vertices))

    total = Rational{Int}(0)
    for dp in _degree_partitions(n_ext, loops, vd)
        degs = sort([d for (d, c) in dp.counts if c > 0])
        isempty(degs) && continue
        cv = Int8[get(dp.counts, d, 0) for d in degs[1]:degs[end]]
        qp = Partition(Int8(n_ext), cv, Int8(degs[1]), Int8(loops))
        s  = TopoState(qp)
        qg21_enumerate!(s) do state
            # Topological filter chain (each reject returns immediately).
            onepi      && !is_one_pi(state)      && return
            nosbridge  && !has_no_sbridge(state) && return
            notadpole  && !has_no_tadpole(state) && return
            onshell    && !has_no_onshell(state) && return
            # qgraf nosn > 0 sets intf(nsl)=1 AND intf(nsb)=1 (f08:2794-2798):
            # nosnail implies both no-self-loop AND no-sbridge.
            if nosnail
                !has_no_snail(state)   && return
                !has_no_sbridge(state) && return
            end
            onevi      && !is_one_vi(state)      && return
            noselfloop && !has_no_selfloop(state) && return
            nodiloop   && !has_no_diloop(state)  && return
            noparallel && !has_no_parallel(state) && return
            labels = compute_qg10_labels(state)
            autos  = enumerate_topology_automorphisms(state)
            g_size = length(autos)
            ps1 = collect(1:n_ext)
            while true
                ext_perm = Symbol[ext_exp[ps1[i]] for i in 1:n_ext]
                qgen_enumerate_assignments(state, labels, ext_perm, dpntro,
                                            exp.conjugate) do st, pmap
                    n_stab = emission_stabilizer(st, labels, autos, ps1, pmap)
                    total += n_stab // g_size
                end
                j = n_ext - 1
                while j >= 1 && ps1[j] >= ps1[j + 1]; j -= 1; end
                j == 0 && break
                k = n_ext
                while ps1[k] <= ps1[j]; k -= 1; end
                ps1[j], ps1[k] = ps1[k], ps1[j]
                reverse!(view(ps1, (j + 1):n_ext))
            end
        end
    end
    @assert denominator(total) == 1 "Burnside sum non-integer ($total)"
    Int(numerator(total))
end

"""
    count_dedup_prefilter(model, in_fields, out_fields; loops) -> Int

Approach (C) — pre-filter ext-leg perms: for each topology, compute
canonical perm representatives (one per orbit under topology auto action
on ext indices) and only enumerate qgen for those.
"""
function count_dedup_prefilter(model, in_fields::Vector{Symbol},
                                 out_fields::Vector{Symbol}; loops::Int)
    n_ext   = length(in_fields) + length(out_fields)
    ext_raw = vcat(in_fields, out_fields)
    exp     = _expand_model_for_diagen(model)
    ext_exp = _expand_external_fields(ext_raw, exp)
    dpntro  = build_dpntro(exp.vertex_rules)
    rules   = feynman_rules(model)
    vd      = Set(length(k) for k in keys(rules.vertices))
    total   = 0
    for dp in _degree_partitions(n_ext, loops, vd)
        degs = sort([d for (d, c) in dp.counts if c > 0])
        isempty(degs) && continue
        mrho = degs[1]
        cv   = Int8[get(dp.counts, d, 0) for d in mrho:degs[end]]
        qp   = Partition(Int8(n_ext), cv, Int8(mrho), Int8(loops))
        s    = TopoState(qp)
        qg21_enumerate!(s) do state
            labels = compute_qg10_labels(state)
            autos  = enumerate_topology_automorphisms(state)
            # Compute canonical ext-perm reps under autos (action on 1..n_ext).
            seen = Set{Vector{Int}}()
            canonical_perms = Vector{Vector{Int}}()
            for ps1 in _permutations(1:n_ext)
                ps1v = collect(ps1)
                ps1v in seen && continue
                # Build orbit under the right-action `(g·ps1)[i] = ps1[g⁻¹[i]]`
                # used elsewhere in this file (`is_emission_canonical`,
                # `_ps1_preserved`, `same_emission_orbit`). Topology autos
                # relabel slots, not physical legs — see the right-action
                # comment block above _ps1_preserved.
                orbit = Set{Vector{Int}}()
                for auto in autos
                    ginv = _inv_perm(auto)
                    new_ps1 = [Int(ps1v[Int(ginv[i])]) for i in 1:n_ext]
                    push!(orbit, new_ps1)
                end
                push!(canonical_perms, minimum(orbit))
                union!(seen, orbit)
            end
            for ps1 in canonical_perms
                ext_perm = Symbol[ext_exp[ps1[i]] for i in 1:n_ext]
                qgen_enumerate_assignments(state, labels, ext_perm, dpntro,
                                            exp.conjugate) do _, _
                    total += 1
                end
            end
        end
    end
    total
end
