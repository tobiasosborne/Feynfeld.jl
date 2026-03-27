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

# Vertex Lorentz structure — dispatch on gauge group + species
# QED: fermion-fermion-photon → γ^μ
function vertex_structure(::U1, ::Fermion, ::Boson, mu::LorentzIndex)
    DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))
end

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
    # Determine species of each field
    f1 = get_field(r.model, fields[1])
    f3 = get_field(r.model, fields[3])
    vertex_structure(r._gauge, species(f1), species(f3), mu)
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
