# Feynfeld.jl — FeynAmpDenominator and propagator types
#
# FAD wraps a product of propagator denominators representing
# 1 / (prop1 * prop2 * ...) in a loop integral.
#
# Three propagator types (mirroring FeynCalc):
#   PropagatorDenominator    — legacy: 1/(q² - m²)
#   StandardPropagator       — modern: 1/((q+p)² + 2q·r - m²)
#   GenericPropagator        — arbitrary: 1/expr
#
# Ref: FeynCalc SharedObjects.m (FeynAmpDenominator, PropagatorDenominator,
#      StandardPropagatorDenominator, GenericPropagatorDenominator)

export FeynAmpDenominator, FAD
export PropagatorDenominator, StandardPropagator, GenericPropagator

"""Abstract supertype for all propagator denominator types."""
abstract type AbstractPropagator end

"""
    PropagatorDenominator(momentum, mass)

Legacy propagator: 1/(q² - m²).
"""
struct PropagatorDenominator <: AbstractPropagator
    momentum::Any   # Momentum or symbolic expression
    mass::Any       # mass (not squared)
end

"""
    StandardPropagator(momentum, eikonal, mass_sq, power, eta_sign)

Modern propagator: 1/((q+p)² + 2*eikonal - m² + i*η)^power.
"""
struct StandardPropagator <: AbstractPropagator
    momentum::Any   # quadratic part (Momentum expression)
    eikonal::Any    # linear/eikonal part (Pair or 0)
    mass_sq::Any    # mass squared (stored with sign: -m² internally)
    power::Int      # propagator power (default 1)
    eta_sign::Int   # +1 (Feynman) or -1
end

StandardPropagator(mom, mass_sq) = StandardPropagator(mom, 0, mass_sq, 1, 1)

"""
    GenericPropagator(expr, power)

Fully general propagator: 1/expr^power.
"""
struct GenericPropagator <: AbstractPropagator
    expr::Any
    power::Int
end

GenericPropagator(expr) = GenericPropagator(expr, 1)

"""
    FeynAmpDenominator(propagators...)

Product of propagator denominators: 1/(prop1 * prop2 * ...).
Represents the denominator of a Feynman amplitude.
"""
struct FeynAmpDenominator <: FeynExpr
    propagators::Vector{AbstractPropagator}
end

FeynAmpDenominator(ps::AbstractPropagator...) = FeynAmpDenominator(collect(ps))

# Equality and hash for propagator types
Base.:(==)(a::PropagatorDenominator, b::PropagatorDenominator) =
    a.momentum == b.momentum && a.mass == b.mass
Base.hash(p::PropagatorDenominator, h::UInt) =
    hash(p.mass, hash(p.momentum, hash(:PD, h)))

Base.:(==)(a::StandardPropagator, b::StandardPropagator) =
    a.momentum == b.momentum && a.eikonal == b.eikonal &&
    a.mass_sq == b.mass_sq && a.power == b.power && a.eta_sign == b.eta_sign
Base.hash(p::StandardPropagator, h::UInt) =
    hash(p.eta_sign, hash(p.power, hash(p.mass_sq, hash(p.eikonal, hash(p.momentum, hash(:SP, h))))))

Base.:(==)(a::GenericPropagator, b::GenericPropagator) =
    a.expr == b.expr && a.power == b.power
Base.hash(p::GenericPropagator, h::UInt) =
    hash(p.power, hash(p.expr, hash(:GP, h)))

Base.:(==)(a::FeynAmpDenominator, b::FeynAmpDenominator) =
    a.propagators == b.propagators
Base.hash(f::FeynAmpDenominator, h::UInt) =
    hash(f.propagators, hash(:FAD, h))

# Multiply FADs = concatenate propagator lists
Base.:*(a::FeynAmpDenominator, b::FeynAmpDenominator) =
    FeynAmpDenominator(vcat(a.propagators, b.propagators))

Base.length(f::FeynAmpDenominator) = length(f.propagators)

"""
    FAD(args...)

Convenience constructor for FeynAmpDenominator.
Each argument is either a Symbol (massless propagator) or a Tuple (momentum, mass).

# Examples
```julia
FAD(:q)              # 1/q²
FAD(:q, :m)          # 1/(q² - m²)
FAD((:q, :m), (:q_minus_p, :m))  # 1/((q²-m²)((q-p)²-m²))
```
"""
function FAD(args...)
    props = AbstractPropagator[]
    for a in args
        if a isa Symbol
            push!(props, PropagatorDenominator(Momentum(a), 0))
        elseif a isa Tuple && length(a) == 2
            push!(props, PropagatorDenominator(Momentum(a[1]), a[2]))
        else
            error("FAD: each arg must be a Symbol or (momentum, mass) tuple")
        end
    end
    FeynAmpDenominator(props)
end
