# Core types: Index, Momentum, Dimension slots
# Design: parametric types + multiple dispatch everywhere.

# ---- Dimension slots (singleton dispatch, unchanged — already Julian) ----
abstract type DimSlot end
struct Dim4   <: DimSlot end
struct DimD   <: DimSlot end
struct DimDm4 <: DimSlot end

dim_trace(::Dim4)   = 4//1
dim_trace(::DimD)   = DIM           # returns DimPoly, not Expr
dim_trace(::DimDm4) = DIM_MINUS_4   # returns DimPoly, not Expr

# BMHV dim_contract: 3x3 dispatch table
dim_contract(::Dim4, ::Dim4)     = Dim4()
dim_contract(::Dim4, ::DimD)     = Dim4()
dim_contract(::DimD, ::Dim4)     = Dim4()
dim_contract(::DimD, ::DimD)     = DimD()
dim_contract(::Dim4, ::DimDm4)   = nothing   # vanishes
dim_contract(::DimDm4, ::Dim4)   = nothing
dim_contract(::DimDm4, ::DimDm4) = DimDm4()
dim_contract(::DimDm4, ::DimD)   = DimDm4()
dim_contract(::DimD, ::DimDm4)   = DimDm4()

# ---- Abstract physics index (enables unified contraction engine) ----
abstract type PhysicsIndex end

# ---- Lorentz index ----
struct LorentzIndex <: PhysicsIndex
    name::Symbol
    dim::DimSlot
end
LorentzIndex(name::Symbol) = LorentzIndex(name, Dim4())

Base.:(==)(a::LorentzIndex, b::LorentzIndex) = a.name == b.name && typeof(a.dim) == typeof(b.dim)
Base.hash(a::LorentzIndex, h::UInt) = hash(a.name, hash(typeof(a.dim), h))
Base.isless(a::LorentzIndex, b::LorentzIndex) = isless(a.name, b.name)
Base.show(io::IO, li::LorentzIndex) = print(io, li.name)

# ---- Momentum ----
struct Momentum
    name::Symbol
    dim::DimSlot
end
Momentum(name::Symbol) = Momentum(name, Dim4())

Base.:(==)(a::Momentum, b::Momentum) = a.name == b.name && typeof(a.dim) == typeof(b.dim)
Base.hash(a::Momentum, h::UInt) = hash(a.name, hash(typeof(a.dim), h))
Base.isless(a::Momentum, b::Momentum) = isless(a.name, b.name)
Base.show(io::IO, m::Momentum) = print(io, m.name)

# ---- MomentumSum (stores full Momentum objects, not Symbols) ----
struct MomentumSum
    terms::Vector{Tuple{Rational{Int}, Momentum}}

    function MomentumSum(terms_raw::Vector{Tuple{Rational{Int}, Momentum}})
        # Combine like terms, sort, drop zeros
        d = Dict{Momentum, Rational{Int}}()
        for (c, m) in terms_raw
            d[m] = get(d, m, Rational{Int}(0)) + c
        end
        terms = sort!([(c, m) for (m, c) in d if !iszero(c)]; by=last)
        isempty(terms) && return nothing  # zero momentum
        length(terms) == 1 && isone(first(terms)[1]) && return first(terms)[2]
        new(terms)
    end
end

Base.show(io::IO, ms::MomentumSum) = join(io, ["$(c)*$(m)" for (c,m) in ms.terms], " + ")

# Arithmetic on momenta
Base.:+(a::Momentum, b::Momentum) = MomentumSum([(1//1, a), (1//1, b)])
Base.:-(a::Momentum, b::Momentum) = MomentumSum([(1//1, a), (-1//1, b)])
Base.:-(m::Momentum) = MomentumSum([(-1//1, m)])
Base.:*(c::Number, m::Momentum) = MomentumSum([(Rational{Int}(c), m)])
Base.:*(m::Momentum, c::Number) = c * m

# Union of things that can index into a Pair
const PairArg = Union{LorentzIndex, Momentum, MomentumSum}
