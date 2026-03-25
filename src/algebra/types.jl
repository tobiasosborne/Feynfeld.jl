# Feynfeld.jl — Algebra layer type definitions
#
# These types form the core type system. Everything above Algebra
# constructs expressions in these types; everything below consumes them.
#
# Phase 0: Foundation types that specialise TensorGR.jl's index machinery
# to the Minkowski metric and QFT-specific algebraic objects.

export LorentzIndex, FourMomentum, MetricTensor
export DiracGamma, DiracChain, SpinorU, SpinorV, Slash
export SUNMatrix, SUNF, SUND, ColourDelta
export PaVe, A0, B0, B1, C0, D0
export Amplitude, FTerm, ScalarProduct, DiracGamma5, FeynExpr

# ── Abstract base ────────────────────────────────────────────────────

"""Abstract supertype for all Feynfeld algebraic expressions."""
abstract type FeynExpr end

# ── Lorentz algebra ──────────────────────────────────────────────────

"""
    LorentzIndex(name, dim=:D)

A Lorentz index μ, ν, ... living in `dim` dimensions.
`dim = :D` for dimensional regularisation, `dim = 4` for 4D.
"""
struct LorentzIndex
    name::Symbol
    dim::Union{Symbol,Int}
end
LorentzIndex(name::Symbol) = LorentzIndex(name, :D)

"""
    FourMomentum(name)

A four-momentum vector p, k, q, ...
"""
struct FourMomentum <: FeynExpr
    name::Symbol
end

"""
    MetricTensor(μ, ν)

The metric tensor g^{μν} (Minkowski signature +−−−).
"""
struct MetricTensor <: FeynExpr
    i::LorentzIndex
    j::LorentzIndex
end

# ── Dirac algebra ────────────────────────────────────────────────────

"""
    DiracGamma(index)

A Dirac gamma matrix γ^μ.
"""
struct DiracGamma <: FeynExpr
    index::LorentzIndex
end

"""
    DiracGamma5()

The chiral matrix γ⁵ = iγ⁰γ¹γ²γ³.
"""
struct DiracGamma5 <: FeynExpr end

"""
    Slash(p)

Feynman slash notation: p̸ = γ^μ p_μ.
"""
struct Slash <: FeynExpr
    momentum::FourMomentum
end

"""
    SpinorU(p, m)

Dirac spinor u(p) for particle with momentum p and mass m.
"""
struct SpinorU <: FeynExpr
    momentum::FourMomentum
    mass::Any
end

"""
    SpinorV(p, m)

Dirac spinor v(p) for antiparticle with momentum p and mass m.
"""
struct SpinorV <: FeynExpr
    momentum::FourMomentum
    mass::Any
end

"""
    DiracChain(elements)

An ordered chain of Dirac-algebra objects (spinors, gammas, slashes)
forming a spinor bilinear like ū(p) γ^μ u(k).
"""
struct DiracChain <: FeynExpr
    elements::Vector{FeynExpr}
end

# ── Colour algebra (SU(N)) ──────────────────────────────────────────

"""
    SUNMatrix(name, a)

A generator T^a of SU(N) in the fundamental representation.
"""
struct SUNMatrix <: FeynExpr
    name::Symbol
    a::Symbol  # colour index
end

"""
    SUNF(a, b, c)

SU(N) structure constants f^{abc}.
"""
struct SUNF <: FeynExpr
    a::Symbol
    b::Symbol
    c::Symbol
end

"""
    SUND(a, b, c)

SU(N) symmetric structure constants d^{abc}.
"""
struct SUND <: FeynExpr
    a::Symbol
    b::Symbol
    c::Symbol
end

"""
    ColourDelta(a, b)

Kronecker delta δ^{ab} in colour space.
"""
struct ColourDelta <: FeynExpr
    a::Symbol
    b::Symbol
end

# ── Passarino-Veltman integral symbols ───────────────────────────────

"""
    PaVe(indices, masses, momenta)

General Passarino-Veltman integral symbol.
Specific cases: A0, B0, B1, C0, D0.
"""
struct PaVe <: FeynExpr
    indices::Vector{Int}
    masses::Vector{Any}
    momenta::Vector{Any}
end

A0(m) = PaVe(Int[], [m], Any[])
B0(p, m1, m2) = PaVe([0], [m1, m2], [p])
B1(p, m1, m2) = PaVe([1], [m1, m2], [p])
C0(p1, p2, p12, m1, m2, m3) = PaVe([0, 0], [m1, m2, m3], [p1, p2, p12])
D0(args...) = PaVe([0, 0, 0], collect(args[4:end]), collect(args[1:3]))

# ── Amplitude ────────────────────────────────────────────────────────

"""
    ScalarProduct(p, q)

Minkowski scalar product p · q = p^μ q_μ.
"""
struct ScalarProduct <: FeynExpr
    p::FourMomentum
    q::FourMomentum
end

"""
    FTerm(coeff, factors)

A single term in an amplitude: a numeric coefficient times
a product of FeynExpr factors.
"""
struct FTerm <: FeynExpr
    coeff::Number
    factors::Vector{FeynExpr}
end

"""
    Amplitude(terms)

A scattering amplitude: a sum of FTerms.
"""
struct Amplitude <: FeynExpr
    terms::Vector{FTerm}
end
