# Layer 1: QFT Model — abstract interface + concrete QED.
#
# Design patterns used:
#   - Abstract interface (MultivariatePolynomials style)
#   - Field species as type parameter (DiracGamma{S} pattern)
#   - Gauge groups as types with compile-time dispatch
#   - Holy traits for orthogonal field properties

# ---- Gauge groups as types ----
abstract type GaugeGroup end
struct U1 <: GaugeGroup end
struct SU{N} <: GaugeGroup
    function SU{N}() where N
        N isa Int && N >= 2 || error("SU(N) requires N ≥ 2, got $N")
        new{N}()
    end
end

# Compile-time Casimirs via dispatch
casimir_fund(::Type{SU{N}}) where N = (N^2 - 1) // (2N)
casimir_adj(::Type{SU{N}}) where N = Rational{Int}(N)
dim_fund(::Type{SU{N}}) where N = N
dim_adj(::Type{SU{N}}) where N = N^2 - 1
Base.show(io::IO, ::U1) = print(io, "U(1)")
Base.show(io::IO, ::SU{N}) where N = print(io, "SU($N)")

# ---- Field species as type parameter ----
abstract type FieldSpecies end
struct Fermion <: FieldSpecies end
struct Boson   <: FieldSpecies end  # spin-1 vector boson
struct Scalar  <: FieldSpecies end

struct Field{S<:FieldSpecies}
    name::Symbol
    mass::Symbol
    charge::Rational{Int}
    self_conjugate::Bool
end

# Convenience constructors
fermion(name, mass; charge=0//1) = Field{Fermion}(name, mass, charge, false)
vector_boson(name, mass; self_conj=true) = Field{Boson}(name, mass, 0//1, self_conj)
scalar(name, mass; charge=0//1) = Field{Scalar}(name, mass, charge, false)

# ---- Holy traits for orthogonal properties ----
struct Massive end
struct Massless end
mass_trait(f::Field) = f.mass == :zero ? Massless() : Massive()

struct Charged end
struct Neutral end
charge_trait(f::Field) = iszero(f.charge) ? Neutral() : Charged()

# Species dispatch
species(::Field{S}) where S = S()

Base.show(io::IO, f::Field{Fermion}) = print(io, "ψ($(f.name))")
Base.show(io::IO, f::Field{Boson}) = print(io, "A($(f.name))")
Base.show(io::IO, f::Field{Scalar}) = print(io, "φ($(f.name))")

# ---- Abstract model interface ----
abstract type AbstractModel end

# Required interface (concrete types MUST implement):
model_name(m::AbstractModel) = error("implement model_name")
model_fields(m::AbstractModel) = error("implement model_fields")
gauge_groups(m::AbstractModel) = error("implement gauge_groups")
model_params(m::AbstractModel) = error("implement model_params")

# Derived interface (free from required methods):
fermion_fields(m::AbstractModel) = [f for f in model_fields(m) if f isa Field{Fermion}]
boson_fields(m::AbstractModel) = [f for f in model_fields(m) if f isa Field{Boson}]
function get_field(m::AbstractModel, name::Symbol)
    idx = findfirst(f -> f.name == name, model_fields(m))
    idx === nothing && error("Field :$name not in model :$(model_name(m))")
    model_fields(m)[idx]
end

# ---- Concrete QED model ----
struct QEDModel <: AbstractModel
    electron::Field{Fermion}
    muon::Field{Fermion}
    photon::Field{Boson}
    params::Dict{Symbol, Any}
end

model_name(::QEDModel) = :QED
model_fields(m::QEDModel) = Field[m.electron, m.muon, m.photon]
gauge_groups(::QEDModel) = GaugeGroup[U1()]
model_params(m::QEDModel) = m.params

function qed_model(; m_e=:m_e, m_mu=:m_mu)
    QEDModel(
        fermion(:e, m_e; charge=-1//1),
        fermion(:mu, m_mu; charge=-1//1),
        vector_boson(:gamma, :zero),
        Dict{Symbol,Any}(:e => :e, :alpha => :(e^2/(4π)))
    )
end

Base.show(io::IO, ::QEDModel) = print(io, "QEDModel(e, μ, γ)")
