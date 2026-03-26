# Feynfeld.jl — Algebra layer type definitions
#
# Core type system. Phase 0 foundation types + Dirac/Colour/PaVe placeholders.
# Lorentz algebra (Pair, Momentum) is in dedicated files.

export LorentzIndex, FourMomentum
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
# Defined in colour_types.jl (Phase 1c)

# ── Passarino-Veltman integral symbols ───────────────────────────────
# Defined in integrals/pave.jl (Phase 1d)

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
