# Layer 2: Feynman rules — callable struct with dispatch on field species.
#
# Design patterns:
#   - Callable struct (rules(:e, :e, :gamma) → vertex factor)
#   - Dispatch on Field{Species} for propagator numerators
#   - Vertex structure encoded via dispatch, not data

# ---- Vertex rule ----
struct VertexRule
    fields::NTuple{3, Symbol}
    coupling::Symbol
end
Base.show(io::IO, v::VertexRule) = print(io, join(v.fields, "-"), " [g=$(v.coupling)]")

# Vertex Lorentz structure — dispatch on gauge group + species + coupling key.
# Returns DiracExpr (may be sum of chains for chiral vertices).

# QED: fermion-fermion-photon → γ^μ
vertex_structure(::GaugeGroup, ::Fermion, ::Boson, ::Val{:e}, mu::LorentzIndex) =
    DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))

# eeZ: neutral current → (g_V - g_A γ5) γ^μ
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Table 10.3
# Note: ordering as (gV - gA γ5)γ^μ is used throughout for internal consistency.
# The γ5-ordering difference vs γ^μ(gV - gA γ5) only affects cross-terms with
# non-chiral channels; the Grozin formula absorbs this convention implicitly.
function vertex_structure(::GaugeGroup, ::Fermion, ::Boson, ::Val{:e_Z}, mu::LorentzIndex)
    gamma_mu = DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))
    g5_gamma_mu = DiracExpr(DiracChain([GA5(), DiracGamma(LISlot(mu))]))
    EW_GV_E_R * gamma_mu - EW_GA_E_R * g5_gamma_mu
end

# eνW: charged current → γ^μ (1-γ5)/2
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eq. (10.63)
# Standard convention: γ^μ P_L (projector RIGHT of gamma).
# P_L γ^μ = γ^μ P_R ≠ γ^μ P_L — ordering matters for s×t interference (Eps term).
vertex_structure(::GaugeGroup, ::Fermion, ::Boson, ::Val{:g_W}, mu::LorentzIndex) =
    DiracExpr(DiracChain([DiracGamma(LISlot(mu)), GA7()]))

# QCD: qqg → γ^μ (same Lorentz structure as QED)
vertex_structure(::GaugeGroup, ::Fermion, ::Boson, ::Val{:g_s}, mu::LorentzIndex) =
    DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))

# Fallback for unknown couplings
vertex_structure(g::GaugeGroup, s1::FieldSpecies, s2::FieldSpecies, ::Val, mu::LorentzIndex) =
    error("No vertex_structure for ($g, $s1, $s2) — implement dispatch")

# ---- Propagator numerators — dispatch on field species ----
# Fermion: i(p-slash + m) / (p² - m²). We return numerator only.
function propagator_num(::Fermion, p::Momentum, mass_val)
    if iszero(mass_val)
        DiracExpr(DiracChain([GS(p)]))
    else
        DiracExpr(DiracChain([GS(p)])) + mass_val * DiracExpr(alg(1))
    end
end

# Vector boson: -ig^{μν}/p² (Feynman gauge). Numerator = -g^{μν}.
function propagator_num(::Boson, mu::LorentzIndex, nu::LorentzIndex)
    -1 * alg(pair(mu, nu))
end

# Scalar: i/(p² - m²). Numerator = 1.
propagator_num(::Scalar, ::Momentum, mass_val) = alg(1)

# ---- FeynmanRules: callable struct ----
struct FeynmanRules
    model::AbstractModel
    vertices::Dict{NTuple{3,Symbol}, VertexRule}
    _gauge::GaugeGroup
end

# Callable: rules(field_names...) → VertexRule
function (r::FeynmanRules)(fields::NTuple{3,Symbol})
    haskey(r.vertices, fields) || error("No vertex for $fields in $(model_name(r.model))")
    r.vertices[fields]
end

# Get the Lorentz structure for a vertex at a given index
function vertex_factor(r::FeynmanRules, fields::NTuple{3,Symbol}, mu::LorentzIndex)
    v = r(fields)
    f1 = get_field(r.model, fields[1])
    f3 = get_field(r.model, fields[3])
    vertex_structure(r._gauge, species(f1), species(f3), Val(v.coupling), mu)
end

function Base.show(io::IO, r::FeynmanRules)
    nv = length(r.vertices)
    print(io, "FeynmanRules($(model_name(r.model)): $(nv) vertices)")
end

# ---- Extract rules from a model ----
function feynman_rules(model::AbstractModel)
    gs = gauge_groups(model)
    length(gs) == 1 || error("Multi-gauge models not yet supported")
    gauge = gs[1]

    vertices = Dict{NTuple{3,Symbol}, VertexRule}()
    for f in fermion_fields(model)
        for b in boson_fields(model)
            key = (f.name, f.name, b.name)
            vertices[key] = VertexRule(key, :e)
        end
    end

    FeynmanRules(model, vertices, gauge)
end
