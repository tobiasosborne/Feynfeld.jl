# Parametric DiracGamma{S} and DiracChain.
# Key change: slot type is a type parameter, enabling dispatch instead of isa cascades.

# ---- Slot types ----
abstract type DiracSlot end

struct LISlot <: DiracSlot
    index::LorentzIndex
end
struct MomSlot <: DiracSlot
    mom::Momentum
end
struct MomSumSlot <: DiracSlot
    mom::MomentumSum
end
struct Gamma5Slot <: DiracSlot end
struct ProjPSlot <: DiracSlot end   # (1+g5)/2
struct ProjMSlot <: DiracSlot end   # (1-g5)/2
struct IdSlot    <: DiracSlot end   # identity

# ---- Parametric DiracGamma ----
struct DiracGamma{S<:DiracSlot}
    slot::S
end

# Convenience constructors
GA(mu::Symbol)  = DiracGamma(LISlot(LorentzIndex(mu)))
GAD(mu::Symbol) = DiracGamma(LISlot(LorentzIndex(mu, DimD())))
GS(p::Symbol)   = DiracGamma(MomSlot(Momentum(p)))
GS(m::Momentum) = DiracGamma(MomSlot(m))
GA5()           = DiracGamma(Gamma5Slot())
GA6()           = DiracGamma(ProjPSlot())
GA7()           = DiracGamma(ProjMSlot())

Base.show(io::IO, g::DiracGamma{LISlot}) = print(io, "γ($(g.slot.index))")
Base.show(io::IO, g::DiracGamma{MomSlot}) = print(io, "γ($(g.slot.mom))")
Base.show(io::IO, g::DiracGamma{Gamma5Slot}) = print(io, "γ5")
Base.show(io::IO, g::DiracGamma{ProjPSlot}) = print(io, "γ+")
Base.show(io::IO, g::DiracGamma{ProjMSlot}) = print(io, "γ-")
Base.show(io::IO, g::DiracGamma{IdSlot}) = print(io, "I")
Base.show(io::IO, g::DiracGamma{MomSumSlot}) = print(io, "γ($(g.slot.mom))")
Base.show(io::IO, g::DiracGamma) = print(io, "γ(?)")

# ---- Spinor kinds as types (not symbols) ----
abstract type SpinorKind end
struct UKind    <: SpinorKind end
struct VKind    <: SpinorKind end
struct UBarKind <: SpinorKind end
struct VBarKind <: SpinorKind end

struct Spinor{K<:SpinorKind}
    momentum::Momentum
    mass::Rational{Int}
end

# Convenience
u(p::Momentum, m=0//1)    = Spinor{UKind}(p, Rational{Int}(m))
v(p::Momentum, m=0//1)    = Spinor{VKind}(p, Rational{Int}(m))
ubar(p::Momentum, m=0//1) = Spinor{UBarKind}(p, Rational{Int}(m))
vbar(p::Momentum, m=0//1) = Spinor{VBarKind}(p, Rational{Int}(m))

Base.show(io::IO, ::Spinor{UKind}) = print(io, "u")
Base.show(io::IO, ::Spinor{VKind}) = print(io, "v")
Base.show(io::IO, ::Spinor{UBarKind}) = print(io, "ū")
Base.show(io::IO, ::Spinor{VBarKind}) = print(io, "v̄")

# Conjugate pairing: dispatch tells us which pairs form spin sums
is_conjugate(::Spinor{UBarKind}, ::Spinor{UKind}) = true
is_conjugate(::Spinor{VBarKind}, ::Spinor{VKind}) = true
is_conjugate(::Spinor, ::Spinor) = false

# ---- DiracChain: ordered product of gamma matrices between spinors ----
const DiracElement = Union{DiracGamma, Spinor}

struct DiracChain
    elements::Vector{DiracElement}
end

# dot() constructor: builds a chain from variadic args, flattening nested chains
function dot(args...)
    elems = DiracElement[]
    for a in args
        if a isa DiracChain
            append!(elems, a.elements)
        else
            push!(elems, a)
        end
    end
    DiracChain(elems)
end

Base.:*(a::DiracChain, b::DiracChain) = DiracChain(vcat(a.elements, b.elements))

Base.:(==)(a::DiracChain, b::DiracChain) = a.elements == b.elements
Base.hash(c::DiracChain, h::UInt) = hash(c.elements, hash(:DiracChain, h))

function Base.show(io::IO, c::DiracChain)
    join(io, c.elements, " . ")
end

# Extract gammas (non-spinor elements) from a chain
gammas(c::DiracChain) = [e for e in c.elements if e isa DiracGamma]

# Extract the Lorentz index from a gamma (dispatch on slot type)
lorentz_index(g::DiracGamma{LISlot}) = g.slot.index
lorentz_index(::DiracGamma) = nothing

# gamma_pair: metric contraction of two gammas → AlgSum (uniform return type).
# Returns zero AlgSum for non-matching pairs (g5, projectors, MomSumSlot).
function _gamma_pair_to_alg(p)
    p isa Number ? (iszero(p) ? AlgSum() : alg(p)) : alg(p)
end
gamma_pair(a::DiracGamma{LISlot}, b::DiracGamma{LISlot}) = _gamma_pair_to_alg(pair(a.slot.index, b.slot.index))
gamma_pair(a::DiracGamma{LISlot}, b::DiracGamma{MomSlot}) = alg(pair(a.slot.index, b.slot.mom))
gamma_pair(a::DiracGamma{MomSlot}, b::DiracGamma{LISlot}) = alg(pair(a.slot.mom, b.slot.index))
gamma_pair(a::DiracGamma{MomSlot}, b::DiracGamma{MomSlot}) = alg(pair(a.slot.mom, b.slot.mom))

# MomSumSlot: handled by expanding before tracing (see dirac_trace.jl).
# gamma_pair does NOT handle MomSumSlot — it falls through to the fallback.
gamma_pair(::DiracGamma, ::DiracGamma) = AlgSum()  # zero: g5, projectors, MomSumSlot
