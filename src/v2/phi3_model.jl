# Layer 1: φ³ scalar field theory model.
#
# Simplest QFT model: single real scalar field with cubic interaction.
# Used for testing diagram generation (pure graph combinatorics, no
# fermion flow or gauge structure).
#
# Ref: Peskin & Schroeder Ch. 4 (scalar field theory)
# Lagrangian: L = ½(∂φ)² - ½m²φ² - (g/3!)φ³

struct Phi3Model <: AbstractModel
    phi::Field{Scalar}
end

"""
    phi3_model(; mass=:m_phi) → Phi3Model

Construct a φ³ scalar field theory model.
"""
function phi3_model(; mass::Symbol=:m_phi)
    phi = Field{Scalar}(:phi, mass, 0//1, true)  # self-conjugate scalar
    Phi3Model(phi)
end

model_name(::Phi3Model) = :phi3
model_fields(m::Phi3Model) = Field[m.phi]
gauge_groups(::Phi3Model) = GaugeGroup[]

"""
    feynman_rules(m::Phi3Model) → FeynmanRules

φ³ has one vertex: [φ, φ, φ] with coupling :g_phi3.
"""
function feynman_rules(m::Phi3Model)
    verts = Dict{Tuple, VertexRule}()
    verts[(:phi, :phi, :phi)] = VertexRule((:phi, :phi, :phi), :g_phi3)
    FeynmanRules(m, verts, U1())
end
