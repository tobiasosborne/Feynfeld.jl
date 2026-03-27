# SU(N) colour algebra types.
# Index types are subtypes of PhysicsIndex (shared with Lorentz).
# All types enforce canonical ordering at construction.

# ---- Colour indices ----
struct AdjointIndex <: PhysicsIndex
    name::Symbol
end
Base.:(==)(a::AdjointIndex, b::AdjointIndex) = a.name == b.name
Base.hash(a::AdjointIndex, h::UInt) = hash(a.name, hash(:AdjIdx, h))
Base.isless(a::AdjointIndex, b::AdjointIndex) = isless(a.name, b.name)
Base.show(io::IO, a::AdjointIndex) = print(io, a.name)

struct FundIndex <: PhysicsIndex
    name::Symbol
end
Base.:(==)(a::FundIndex, b::FundIndex) = a.name == b.name
Base.hash(a::FundIndex, h::UInt) = hash(a.name, hash(:FundIdx, h))
Base.isless(a::FundIndex, b::FundIndex) = isless(a.name, b.name)
Base.show(io::IO, a::FundIndex) = print(io, a.name)

# ---- Generator T^a (lives inside ColourChain, not an AlgFactor) ----
struct SUNT
    adj::AdjointIndex
end
Base.show(io::IO, t::SUNT) = print(io, "T($(t.adj))")

# ---- Adjoint Kronecker delta δ^{ab} (canonical: a ≤ b) ----
struct SUNDelta
    a::AdjointIndex
    b::AdjointIndex
    function SUNDelta(a::AdjointIndex, b::AdjointIndex)
        isless(a, b) ? new(a, b) : new(b, a)
    end
end
Base.:(==)(x::SUNDelta, y::SUNDelta) = x.a == y.a && x.b == y.b
Base.hash(d::SUNDelta, h::UInt) = hash(d.a, hash(d.b, hash(:SUNDelta, h)))
Base.isless(x::SUNDelta, y::SUNDelta) = x.a == y.a ? isless(x.b, y.b) : isless(x.a, y.a)
Base.show(io::IO, d::SUNDelta) = print(io, "δ($(d.a),$(d.b))")

# ---- Fundamental Kronecker delta δ_{ij} (canonical: i ≤ j) ----
struct FundDelta
    i::FundIndex
    j::FundIndex
    function FundDelta(i::FundIndex, j::FundIndex)
        isless(i, j) ? new(i, j) : new(j, i)
    end
end
Base.:(==)(x::FundDelta, y::FundDelta) = x.i == y.i && x.j == y.j
Base.hash(d::FundDelta, h::UInt) = hash(d.i, hash(d.j, hash(:FundDelta, h)))
Base.isless(x::FundDelta, y::FundDelta) = x.i == y.i ? isless(x.j, y.j) : isless(x.i, y.i)
Base.show(io::IO, d::FundDelta) = print(io, "δ_f($(d.i),$(d.j))")

# ---- Antisymmetric structure constant f^{abc} (sorted, with sign) ----
struct SUNF
    a::AdjointIndex
    b::AdjointIndex
    c::AdjointIndex
    sign::Int  # +1 or -1 from antisymmetric sort
    function SUNF(a::AdjointIndex, b::AdjointIndex, c::AdjointIndex)
        indices = [a, b, c]
        sorted = sort(indices)
        sign = _parity([findfirst(==(s), indices) for s in sorted])
        new(sorted[1], sorted[2], sorted[3], sign)
    end
end
Base.:(==)(x::SUNF, y::SUNF) = x.a == y.a && x.b == y.b && x.c == y.c && x.sign == y.sign
Base.hash(f::SUNF, h::UInt) = hash(f.a, hash(f.b, hash(f.c, hash(f.sign, hash(:SUNF, h)))))
function Base.isless(x::SUNF, y::SUNF)
    x.a != y.a && return isless(x.a, y.a)
    x.b != y.b && return isless(x.b, y.b)
    x.c != y.c && return isless(x.c, y.c)
    return x.sign < y.sign
end
Base.show(io::IO, f::SUNF) = print(io, f.sign < 0 ? "-f" : "f", "($(f.a),$(f.b),$(f.c))")

# ---- Symmetric structure constant d^{abc} (sorted) ----
struct SUND
    a::AdjointIndex
    b::AdjointIndex
    c::AdjointIndex
    function SUND(a::AdjointIndex, b::AdjointIndex, c::AdjointIndex)
        sorted = sort([a, b, c])
        new(sorted[1], sorted[2], sorted[3])
    end
end
Base.:(==)(x::SUND, y::SUND) = x.a == y.a && x.b == y.b && x.c == y.c
Base.hash(d::SUND, h::UInt) = hash(d.a, hash(d.b, hash(d.c, hash(:SUND, h))))
function Base.isless(x::SUND, y::SUND)
    x.a != y.a && return isless(x.a, y.a)
    x.b != y.b && return isless(x.b, y.b)
    return isless(x.c, y.c)
end
Base.show(io::IO, d::SUND) = print(io, "d($(d.a),$(d.b),$(d.c))")

# ---- ColourChain: ordered product of SUNT generators ----
struct ColourChain
    elements::Vector{SUNT}
end
Base.length(c::ColourChain) = length(c.elements)
Base.show(io::IO, c::ColourChain) = join(io, c.elements, ".")

# ---- Permutation parity ----
function _parity(perm::Vector{Int})
    n = length(perm)
    visited = falses(n)
    sign = 1
    for i in 1:n
        visited[i] && continue
        j, len = i, 0
        while !visited[j]
            visited[j] = true
            j = perm[j]
            len += 1
        end
        iseven(len) && (sign = -sign)
    end
    sign
end

# ---- Fresh index counter (for dummy indices in traces) ----
const _COLOUR_DUMMY_COUNTER = Ref(0)
function _fresh_adj()
    _COLOUR_DUMMY_COUNTER[] += 1
    AdjointIndex(Symbol("_c$(_COLOUR_DUMMY_COUNTER[])"))
end
