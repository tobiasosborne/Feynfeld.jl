#  Phase 18a-3: per-edge propagator factors.
#
#  For each internal edge of an emitted topology, look up the field at the
#  edge (via pmap), get its species from the model, and emit a Propagator
#  bundle holding the propagator numerator + denominator AlgSum factors.
#
#  Tree-level scope:
#    Boson  : num = alg(1)  (metric is implicit via shared Lorentz-index
#                            convention — matches src/v2/amplitude.jl)
#    Fermion: num = DiracExpr(p̸ + m) via propagator_num
#    Scalar : num = alg(1) via propagator_num
#  Denominator: AlgSum representing (p² − m²).  For massless fields this is
#  just `alg(pair(mom, mom))`; massive support arrives in 18b.

"""
    Propagator

One per internal edge of an emitted diagram. Carries everything Phase
18a-7 needs to assemble the amplitude for that edge:

- `edge_id`: matches `compute_amap` (= ege[v_lo, v_hi] + parallel_offset).
- `v_lo, v_hi, parallel_idx`: edge endpoints (mirror `InternalEdge`).
- `field`: field name from `pmap` at the v_lo half-edge (the v_hi half
  carries its conjugate; either side gives the same SPECIES).
- `mom`: edge momentum from `route_momenta` (qgraf "all incoming"
  convention — caller-side handles physical sign).
- `num`: propagator numerator. Boson → `alg(1)`; fermion →
  `DiracExpr(p̸+m)`; scalar → `alg(1)`.
- `denom`: AlgSum representing `(p² − m²)`.  Massless: `alg(pair(p, p))`.
"""
struct Propagator
    edge_id::Int
    v_lo::Int
    v_hi::Int
    parallel_idx::Int
    field::Symbol
    mom::Union{Nothing, Momentum, MomentumSum}
    num::Union{AlgSum, DiracExpr}
    denom::AlgSum
end

"""
    build_propagators(state, labels, pmap, edge_mom, model) -> Vector{Propagator}

One Propagator per internal edge in `edge_mom.internal`. Field species
drives numerator dispatch; mass drives denominator construction.

Tree-level only (massless propagators or rational placeholder masses).
PaVe machinery and explicit `-g^{μν}` metric for boson edges are
deferred to Phase 18b.
"""
function build_propagators(state::TopoState, labels,
                            pmap::AbstractMatrix{Symbol},
                            edge_mom::EdgeMomenta,
                            model::AbstractModel)
    propagators = Propagator[]
    amap = compute_amap(state, labels)
    for (k, ie) in enumerate(edge_mom.internal)
        # Find the slot at v_lo connecting to v_hi at this parallel index.
        # vmap rows list neighbours in order; the parallel_idx-th matching
        # slot is the one for this edge.
        slot = _find_slot(labels, ie.v_lo, ie.v_hi, ie.parallel_idx)
        field = pmap[ie.v_lo, slot]
        f = get_field(model, field)
        mass_val = f.mass == :zero ? 0//1 : 1//1   # placeholder for symbolic mass
        num   = _propagator_numerator(species(f), ie.momentum, mass_val)
        denom = _propagator_denominator(ie.momentum, mass_val)
        push!(propagators,
              Propagator(amap[ie.v_lo, slot], ie.v_lo, ie.v_hi,
                         ie.parallel_idx, field, ie.momentum, num, denom))
    end
    propagators
end

# Slot at v_lo whose vmap entry equals v_hi, picking the parallel_idx-th
# such match in slot order.
function _find_slot(labels, v_lo::Int, v_hi::Int, parallel_idx::Int)
    seen = 0
    @inbounds for s in 1:size(labels.vmap, 2)
        Int(labels.vmap[v_lo, s]) == v_hi || continue
        seen += 1
        seen == parallel_idx && return s
    end
    error("build_propagators: no slot at $v_lo matches v_hi=$v_hi at parallel_idx=$parallel_idx")
end

# Numerator: index-sharing convention for bosons (matches amplitude.jl);
# explicit (p̸ + m) for fermions; identity for scalars (propagator_num
# would return the same alg(1) so we inline rather than dispatch through
# the Momentum-only overload).
_propagator_numerator(::Boson,  _mom, _mass) = alg(1)
_propagator_numerator(::Scalar, _mom, _mass) = alg(1)
function _propagator_numerator(::Fermion, mom::Momentum, mass)
    propagator_num(Fermion(), mom, mass)
end
# Composite momentum: γ^μ is linear in its argument, so for p = Σ cᵢ pᵢ
# the numerator (p̸ + m) expands to Σ cᵢ p̸ᵢ + m·I.
# Phase 18b-2 (feynfeld-h3pb) — required for internal fermion propagators
# in Compton/Bhabha/qq̄→gg.
function _propagator_numerator(::Fermion, mom::MomentumSum, mass)
    de = DiracExpr()
    for (c, p) in mom.terms
        de = de + c * DiracExpr(DiracChain([GS(p)]))
    end
    iszero(mass) ? de : de + mass * DiracExpr(alg(1))
end

# Denominator (p² − m²). For massless: just pair(mom, mom).
function _propagator_denominator(mom, mass_val)
    sp = alg(pair(mom, mom))
    iszero(mass_val) && return sp
    sp - mass_val * alg(1)
end
