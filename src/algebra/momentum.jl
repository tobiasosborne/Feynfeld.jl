# Feynfeld.jl — Momentum type for Lorentz algebra
#
# Momentum wraps a name (Symbol) and dimension slot (DimSlot).
# This is the QFT four-momentum — a slot argument for Pair, not a FeynExpr.
#
# Ported from FeynCalc SharedObjects.m Momentum[p, dim].

export Momentum, PairArg

"""
    Momentum(name::Symbol, dim::DimSlot=Dim4())

A four-momentum in the given dimension slot.
Default is 4-dimensional, matching FeynCalc convention `Momentum[p]`.

# Examples
```julia
Momentum(:p)           # 4D momentum
Momentum(:p, DimD())   # D-dimensional momentum
Momentum(:p, DimDm4()) # evanescent (D-4)-dimensional
```
"""
struct Momentum
    name::Symbol
    dim::DimSlot
end

Momentum(name::Symbol) = Momentum(name, Dim4())
Momentum(name::Symbol, d::Int) = Momentum(name, to_dim(d))
Momentum(name::Symbol, d::Symbol) = Momentum(name, to_dim(d))

function Base.isless(a::Momentum, b::Momentum)
    a.name < b.name || (a.name == b.name && _dim_order(a.dim) < _dim_order(b.dim))
end

Base.hash(m::Momentum, h::UInt) = hash(m.dim, hash(m.name, hash(:Momentum, h)))

"""
    PairArg

Union of valid slot arguments for Pair: LorentzIndex or Momentum.
MomentumSum will be added in a later phase for ExpandScalarProduct support.
"""
const PairArg = Union{LorentzIndex, Momentum}
