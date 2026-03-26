# Feynfeld.jl — Algebra layer type definitions
#
# Core type system. Phase 0 foundation types + Dirac/Colour/PaVe placeholders.
# Lorentz algebra (Pair, Momentum) is in dedicated files.

export LorentzIndex, FourMomentum
export SUNMatrix, SUNF, SUND, ColourDelta
export PaVe, A0, B0, B1, C0, D0
export Amplitude, FTerm, FeynExpr

# ── Abstract base ────────────────────────────────────────────────────

"""Abstract supertype for all Feynfeld algebraic expressions."""
abstract type FeynExpr end

# ── Lorentz index ────────────────────────────────────────────────────

"""
    LorentzIndex(name, dim=Dim4())

A Lorentz index μ, ν, ... living in `dim` dimensions.
Default is 4-dimensional, matching FeynCalc convention.
For dimensional regularisation, use `LorentzIndex(:μ, DimD())`.
"""
struct LorentzIndex
    name::Symbol
    dim::DimSlot
end

LorentzIndex(name::Symbol) = LorentzIndex(name, Dim4())
LorentzIndex(name::Symbol, d::Int) = LorentzIndex(name, to_dim(d))
LorentzIndex(name::Symbol, d::Symbol) = LorentzIndex(name, to_dim(d))

# DimSlot ordering for canonical comparisons
_dim_order(::Dim4) = 1
_dim_order(::DimD) = 2
_dim_order(::DimDm4) = 3

function Base.isless(a::LorentzIndex, b::LorentzIndex)
    a.name < b.name || (a.name == b.name && _dim_order(a.dim) < _dim_order(b.dim))
end

Base.hash(li::LorentzIndex, h::UInt) = hash(li.dim, hash(li.name, hash(:LorentzIndex, h)))

"""
    FourMomentum(name)

Legacy four-momentum type from Phase 0 scaffold.
For new code, prefer `Momentum(name, dim)`.
"""
struct FourMomentum <: FeynExpr
    name::Symbol
end

# ── Dirac algebra ────────────────────────────────────────────────────
# Defined in dirac_types.jl (Phase 1b)

# ── Colour algebra (SU(N)) ──────────────────────────────────────────

"""
    SUNMatrix(name, a)

A generator T^a of SU(N) in the fundamental representation.
"""
struct SUNMatrix <: FeynExpr
    name::Symbol
    a::Symbol
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

# ── Expression tree (Phase 0 placeholder) ────────────────────────────

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
