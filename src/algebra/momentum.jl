# Feynfeld.jl — Momentum types for Lorentz algebra
#
# Momentum wraps a name (Symbol) and dimension slot (DimSlot).
# MomentumSum represents linear combinations of momenta.
# These are slot arguments for Pair, not FeynExpr subtypes.
#
# Ported from FeynCalc SharedObjects.m Momentum[p, dim].

export Momentum, MomentumSum, PairArg

# ── Momentum ─────────────────────────────────────────────────────────

"""
    Momentum(name::Symbol, dim::DimSlot=Dim4())

A four-momentum in the given dimension slot.
Default is 4-dimensional, matching FeynCalc convention `Momentum[p]`.
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

# ── MomentumSum ──────────────────────────────────────────────────────

"""
    MomentumSum(terms, dim)

A linear combination of momenta: sum_i c_i * p_i, all in one dim slot.
Invariant: terms sorted by symbol, no duplicate symbols, no zero coefficients.

Construct via arithmetic: `Momentum(:p) + Momentum(:q)`.
"""
struct MomentumSum
    terms::Vector{Tuple{Rational{Int},Symbol}}
    dim::DimSlot
end

"""Normalize a list of (coeff, name) terms into a MomentumSum, Momentum, or nothing."""
function _normalize_mom(raw::Vector{Tuple{Rational{Int},Symbol}}, dim::DimSlot)
    merged = Dict{Symbol,Rational{Int}}()
    for (c, s) in raw
        merged[s] = get(merged, s, 0 // 1) + c
    end
    filtered = sort!([(c, s) for (s, c) in merged if !iszero(c)]; by=last)
    isempty(filtered) && return nothing
    length(filtered) == 1 && isone(first(filtered)[1]) && return Momentum(first(filtered)[2], dim)
    MomentumSum(filtered, dim)
end

function Base.isless(a::MomentumSum, b::MomentumSum)
    length(a.terms) != length(b.terms) && return length(a.terms) < length(b.terms)
    for (ta, tb) in zip(a.terms, b.terms)
        ta[2] != tb[2] && return ta[2] < tb[2]
        ta[1] != tb[1] && return ta[1] < tb[1]
    end
    false
end

# Cross-type ordering for Pair canonical form
Base.isless(::Momentum, ::MomentumSum) = true
Base.isless(::MomentumSum, ::Momentum) = false

Base.hash(ms::MomentumSum, h::UInt) = hash(ms.dim, hash(ms.terms, hash(:MomentumSum, h)))

# ── Arithmetic ───────────────────────────────────────────────────────

function _check_dim(a, b)
    da, db = a.dim, b.dim
    da === db || error("Cannot combine momenta in different dimensions: $da, $db")
    da
end

function Base.:+(a::Momentum, b::Momentum)
    d = _check_dim(a, b)
    _normalize_mom([(1 // 1, a.name), (1 // 1, b.name)], d)
end

function Base.:-(a::Momentum, b::Momentum)
    d = _check_dim(a, b)
    _normalize_mom([(1 // 1, a.name), (-1 // 1, b.name)], d)
end

Base.:-(a::Momentum) = MomentumSum([(-1 // 1, a.name)], a.dim)

function Base.:*(n::Union{Integer,Rational}, m::Momentum)
    iszero(n) && return nothing
    r = Rational{Int}(n)
    isone(r) && return m
    MomentumSum([(r, m.name)], m.dim)
end
Base.:*(m::Momentum, n::Union{Integer,Rational}) = n * m

function Base.:+(a::MomentumSum, b::Momentum)
    d = _check_dim(a, b)
    _normalize_mom(vcat(a.terms, [(1 // 1, b.name)]), d)
end
Base.:+(a::Momentum, b::MomentumSum) = b + a

function Base.:+(a::MomentumSum, b::MomentumSum)
    d = _check_dim(a, b)
    _normalize_mom(vcat(a.terms, b.terms), d)
end

function Base.:-(a::MomentumSum, b::Momentum)
    d = _check_dim(a, b)
    _normalize_mom(vcat(a.terms, [(-1 // 1, b.name)]), d)
end

function Base.:-(a::Momentum, b::MomentumSum)
    d = _check_dim(a, b)
    _normalize_mom(vcat([(1 // 1, a.name)], [(-c, s) for (c, s) in b.terms]), d)
end

Base.:-(a::MomentumSum) = MomentumSum([(-c, s) for (c, s) in a.terms], a.dim)

function Base.:*(n::Union{Integer,Rational}, ms::MomentumSum)
    iszero(n) && return nothing
    r = Rational{Int}(n)
    MomentumSum([(r * c, s) for (c, s) in ms.terms], ms.dim)
end
Base.:*(ms::MomentumSum, n::Union{Integer,Rational}) = n * ms

# ── PairArg ──────────────────────────────────────────────────────────

"""Union of valid slot arguments for Pair."""
const PairArg = Union{LorentzIndex, Momentum, MomentumSum}
