# Feynfeld.jl — SU(N) colour algebra type system
#
# Two index spaces: adjoint (SUNIndex, 1..N²-1) and fundamental (SUNFIndex, 1..N).
# Core objects: SUNT (generator), SUNTF (explicit-index chain), SUNF/SUND (structure
# constants), SUNDelta/SUNFDelta (Kronecker deltas).
#
# Ref: FeynCalc SharedObjects.m (SUNT, SUNTF, SUNF, SUND, SUNDelta, SUNFDelta)

export SUNIndex, SUNFIndex
export SUNT, SUNTF, SUNF, SUND, SUNDelta, SUNFDelta
export ColourChain, SUNN, CA, CF

# ── Index types ──────────────────────────────────────────────────────

"""Adjoint colour index (range 1..N²-1)."""
struct SUNIndex
    name::Symbol
end

"""Fundamental colour index (range 1..N)."""
struct SUNFIndex
    name::Symbol
end

Base.hash(i::SUNIndex, h::UInt) = hash(i.name, hash(:SUNIndex, h))
Base.hash(i::SUNFIndex, h::UInt) = hash(i.name, hash(:SUNFIndex, h))

# ── Generators and structure constants ───────────────────────────────

"""
    SUNT(a::SUNIndex)

SU(N) generator T^a in the fundamental representation.
Implicit fundamental indices (non-commutative, lives in ColourChain).
"""
struct SUNT <: FeynExpr
    a::SUNIndex
end

"""
    SUNTF(adj, i, j)

Product of generators (T^{a1}...T^{an})_{ij} with explicit fundamental indices.
Empty `adj` vector = SUNFDelta(i,j).
"""
struct SUNTF <: FeynExpr
    adj::Vector{SUNIndex}
    i::SUNFIndex
    j::SUNFIndex
end

"""
    SUNF(a, b, c)

Totally antisymmetric SU(N) structure constants f^{abc}.
Canonical ordering enforced at construction (with sign).
"""
struct SUNF <: FeynExpr
    a::SUNIndex
    b::SUNIndex
    c::SUNIndex
    sign::Int  # +1 or -1 from antisymmetric sorting
    function SUNF(a::SUNIndex, b::SUNIndex, c::SUNIndex)
        args = [a, b, c]
        perm = sortperm(args; by=x -> x.name)
        sgn = _parity(perm)
        new(args[perm[1]], args[perm[2]], args[perm[3]], sgn)
    end
end

"""
    SUND(a, b, c)

Totally symmetric SU(N) structure constants d^{abc}.
Canonically ordered (sorted by name).
"""
struct SUND <: FeynExpr
    a::SUNIndex
    b::SUNIndex
    c::SUNIndex
    function SUND(a::SUNIndex, b::SUNIndex, c::SUNIndex)
        args = sort([a, b, c]; by=x -> x.name)
        new(args[1], args[2], args[3])
    end
end

"""Kronecker delta δ^{ab} in adjoint colour space."""
struct SUNDelta <: FeynExpr
    a::SUNIndex
    b::SUNIndex
    function SUNDelta(a::SUNIndex, b::SUNIndex)
        a.name <= b.name ? new(a, b) : new(b, a)
    end
end

"""Kronecker delta δ_{ij} in fundamental colour space."""
struct SUNFDelta <: FeynExpr
    i::SUNFIndex
    j::SUNFIndex
    function SUNFDelta(i::SUNFIndex, j::SUNFIndex)
        i.name <= j.name ? new(i, j) : new(j, i)
    end
end

# ── ColourChain ──────────────────────────────────────────────────────

"""Non-commutative chain of SUNT generators (colour analog of DiracChain)."""
struct ColourChain <: FeynExpr
    elements::Vector{SUNT}
end

Base.length(c::ColourChain) = length(c.elements)

# ── Casimir invariants ───────────────────────────────────────────────

"""Symbolic number of colours N."""
const SUNN = :N

"""Adjoint Casimir CA = N."""
const CA = :CA

"""Fundamental Casimir CF = (N²-1)/(2N)."""
const CF = :CF

# ── Helpers ──────────────────────────────────────────────────────────

"""Parity of a permutation: +1 for even, -1 for odd."""
function _parity(perm::Vector{Int})
    n = length(perm)
    visited = falses(n)
    sgn = 1
    for i in 1:n
        visited[i] && continue
        visited[i] = true
        j = perm[i]
        cycle_len = 1
        while j != i
            visited[j] = true
            j = perm[j]
            cycle_len += 1
        end
        iseven(cycle_len) && (sgn *= -1)
    end
    sgn
end

Base.isless(a::SUNIndex, b::SUNIndex) = a.name < b.name
Base.isless(a::SUNFIndex, b::SUNFIndex) = a.name < b.name
