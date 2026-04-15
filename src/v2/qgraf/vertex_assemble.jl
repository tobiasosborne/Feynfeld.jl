#  Phase 18a-4: per-vertex Lorentz factors.
#
#  For each internal vertex of an emitted topology, build the Lorentz
#  factor by allocating a SHARED Lorentz index per attached boson edge
#  (matches src/v2/amplitude.jl convention — both endpoints carry the
#  same index name so contraction happens automatically downstream).
#
#  Index naming: per edge id `e`, index = LorentzIndex(Symbol(:mu_l_, e), DimD()).
#  Edge ids come from compute_amap so two endpoints of one edge agree.
#
#  Tree-level scope: 3-vertices with one boson (QED/QCD/EW chiral
#  vertices) or all-scalars (φ³).  4-vertex (gggg) errors deliberately;
#  Phase 18b will add it.

"""
    build_vertices(state, labels, pmap, edge_mom, model) -> Dict{Int, DiracExpr}

For each internal vertex `v` (rhop1..n) build the Lorentz factor as
`DiracExpr`. For QED/QCD-style 3-vertices with one boson leg, returns
`γ^μ` where `μ = LorentzIndex(Symbol(:mu_l_, edge_id), DimD())` and
`edge_id` is the global half-edge id of the boson edge at `v`.

For all-scalar vertices (φ³): returns `DiracExpr(alg(1))`.

For 4-vertices: errors with a Phase-18b deferral message.
"""
function build_vertices(state::TopoState, labels,
                         pmap::AbstractMatrix{Symbol},
                         edge_mom::EdgeMomenta,
                         model::AbstractModel)
    rules = Main.FeynfeldX.feynman_rules(model)
    amap  = compute_amap(state, labels)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)

    out = Dict{Int, DiracExpr}()
    for v in rhop1:n
        vdeg_v = Int(state.vdeg[v])
        fields = Symbol[pmap[v, s] for s in 1:vdeg_v]
        out[v] = _vertex_factor_at(v, fields, vdeg_v, labels, amap, model, rules)
    end
    out
end

function _vertex_factor_at(v::Int, fields::Vector{Symbol}, vdeg_v::Int,
                            labels, amap, model, rules)
    if vdeg_v != 3
        error("build_vertices: vertex $v has degree $vdeg_v; only 3-vertices supported (Phase 18a, 4-vertex deferred to 18b)")
    end

    # pmap fields may carry _bar suffixes (e.g. :e_bar) that aren't in
    # the unexpanded model — canonicalise before querying species.
    species_at = [_field_species(model, fields[s]) for s in 1:vdeg_v]
    boson_slots = Int[s for s in 1:vdeg_v if species_at[s] isa Boson]

    # All-scalar 3-vertex (φ³): no Lorentz structure.
    if isempty(boson_slots) && all(sp -> sp isa Scalar, species_at)
        return DiracExpr(alg(1))
    end

    # QED/QCD/EW chiral 3-vertex: 2 fermions + 1 boson.
    if length(boson_slots) == 1
        slot_b   = boson_slots[1]
        edge_id  = Int(amap[v, slot_b])
        mu       = Main.FeynfeldX.LorentzIndex(Symbol(:mu_l_, edge_id), Main.FeynfeldX.DimD())
        key      = _vertex_rule_key(fields, model)
        haskey(rules.vertices, key) ||
            error("build_vertices: no rule for $key (canonicalised from $fields at v=$v)")
        return Main.FeynfeldX.vertex_factor(rules, key, mu)
    end

    error("build_vertices: 3-vertex at v=$v has $(length(boson_slots)) bosons (only 0 or 1 supported)")
end

# Canonicalize pmap fields (which may carry _bar suffixes for fermion
# antiparticles) to the rule-dict key form: strip _bar from each fermion
# antiparticle name, then arrange as (fermion, fermion, boson) per the
# FeynmanRules.vertices key convention (rules.jl:113-118).
function _vertex_rule_key(fields::Vector{Symbol}, model)
    canon = Symbol[_canonical_field_name(f) for f in fields]
    sort!(canon; by = f -> _field_species(model, f) isa Boson ? 1 : 0)
    Tuple(canon)
end

function _canonical_field_name(f::Symbol)
    s = String(f)
    endswith(s, "_bar") ? Symbol(s[1:end-4]) : f
end

# Look up species for a pmap field; canonicalises _bar suffixes since
# the model dict only stores particles, not antiparticles.
function _field_species(model, name::Symbol)
    Main.FeynfeldX.species(Main.FeynfeldX.get_field(model, _canonical_field_name(name)))
end

# ── Phase 18a-5: per-external spinor / polarisation factors ─────────────

"""
    ExternalFactor

One per external leg. Carries the spinor / polarisation factor and the
in/out + antiparticle metadata Phase 18a-6 (fermion-line traversal)
needs to chain spinors into Dirac chains.

For boson externals: `spinor === nothing` and `position === nothing`.
Polarisation handling is deferred to a future spin-sum stage; if the
boson is internal-only (e.g. tree QED ee→μμ) we never need it.
"""
struct ExternalFactor
    leg_idx::Int
    field::Symbol
    momentum::Momentum
    incoming::Bool
    antiparticle::Bool
    spinor::Union{Nothing, Main.FeynfeldX.Spinor}
    position::Union{Nothing, Symbol}    # :left | :right for fermion
end

"""
    build_externals(state, pmap, physical_moms, n_inco, model) -> Vector{ExternalFactor}

For each external leg i (1..n_ext), construct ExternalFactor:
  - field      = pmap[i, 1]
  - momentum   = physical_moms[i] (PHYSICAL — not qgraf "all incoming")
  - incoming   = i ≤ n_inco
  - antiparticle = field name ends with `_bar`
  - spinor / position via the standard u/v/ubar/vbar dispatch

Mass: read from the model field (e.g. `:zero` or `:m_e`); for `:zero`
spinors get mass=0, otherwise mass=1//1 placeholder (matches
amplitude.jl:115 convention; symbolic masses arrive in 18b).
"""
function build_externals(state::TopoState,
                          pmap::AbstractMatrix{Symbol},
                          physical_moms::Vector{Momentum},
                          n_inco::Int,
                          model::AbstractModel)
    n_ext = Int(state.n_ext)
    length(physical_moms) == n_ext ||
        error("build_externals: physical_moms length $(length(physical_moms)) ≠ n_ext $n_ext")

    out = ExternalFactor[]
    for i in 1:n_ext
        field    = pmap[i, 1]
        mom      = physical_moms[i]
        incoming = i <= n_inco
        anti     = _is_antiparticle_field(field)
        species  = _field_species(model, field)
        spin, pos = _spinor_dispatch(species, incoming, anti, mom,
                                       _ext_mass(model, field))
        push!(out, ExternalFactor(i, field, mom, incoming, anti, spin, pos))
    end
    out
end

_is_antiparticle_field(f::Symbol) = endswith(String(f), "_bar")

function _ext_mass(model, name::Symbol)
    f = Main.FeynfeldX.get_field(model, _canonical_field_name(name))
    f.mass == :zero ? 0//1 : 1//1
end

# Mirror of src/v2/amplitude.jl:167-179 dispatch table for fermions.
function _spinor_dispatch(::Fermion, incoming::Bool, anti::Bool,
                            p::Momentum, m::Rational{Int})
    if      incoming &&  !anti;  (Main.FeynfeldX.u(p, m),    :right)
    elseif  incoming &&   anti;  (Main.FeynfeldX.vbar(p, m), :left)
    elseif !incoming &&  !anti;  (Main.FeynfeldX.ubar(p, m), :left)
    else                          (Main.FeynfeldX.v(p, m),    :right)
    end
end

# Boson externals: polarisation handled at spin-sum stage (deferred).
_spinor_dispatch(::Union{Boson, Scalar}, _, _, _, _) = (nothing, nothing)
