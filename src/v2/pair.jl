# Parametric Pair{A,B}: the universal Lorentz bilinear.
# NOT exported (conflicts with Base.Pair). Users use SP/FV/MT/SPD/FVD/MTD.
#
# Pair(LI,LI) = metric, Pair(LI,M) = four-vector, Pair(M,M) = scalar product.
# Constructor always succeeds (type-stable). BMHV check lives in pair() factory.

struct Pair{A<:PairArg, B<:PairArg}
    a::A
    b::B
    function Pair(a::A, b::B) where {A<:PairArg, B<:PairArg}
        # Canonical ordering only — no BMHV check here (that's pair()'s job)
        if _pair_order(a) > _pair_order(b)
            return new{B, A}(b, a)
        elseif _pair_order(a) == _pair_order(b) && _pair_less(b, a)
            return new{B, A}(b, a)
        end
        new{A, B}(a, b)
    end
end

# Ordering helpers
_pair_order(::LorentzIndex) = 1
_pair_order(::Momentum) = 2
_pair_order(::MomentumSum) = 3

_pair_less(a::LorentzIndex, b::LorentzIndex) = isless(a, b)
_pair_less(a::Momentum, b::Momentum) = isless(a, b)
_pair_less(::PairArg, ::PairArg) = false

# Type aliases — the three physical meanings
const MetricTensor  = Pair{LorentzIndex, LorentzIndex}
const FourVector    = Union{Pair{LorentzIndex, Momentum}, Pair{Momentum, LorentzIndex}}
const ScalarProduct = Pair{Momentum, Momentum}

# Public factory. Always returns a Pair (type-stable; concrete by dispatch).
# BMHV vanishing for mixed-dim LI×LI contractions is the contraction engine's
# concern — see contract.jl::_do_contraction(MT,MT) and
# eps_contract.jl::_det4x4_pairs, both of which call dim_contract themselves.
pair(a::PairArg, b::PairArg) = Pair(a, b)

# ---- Equality, hashing, ordering (structural, not repr-based) ----
Base.:(==)(a::Pair, b::Pair) = a.a == b.a && a.b == b.b
Base.hash(p::Pair, h::UInt) = hash(p.a, hash(p.b, hash(:Pair, h)))

# Structural isless: order by type tag, then by content
function Base.isless(a::Pair, b::Pair)
    oa, ob = _pair_order(a.a), _pair_order(b.a)
    oa != ob && return oa < ob
    if a.a != b.a
        return isless(a.a, b.a)
    end
    oa2, ob2 = _pair_order(a.b), _pair_order(b.b)
    oa2 != ob2 && return oa2 < ob2
    return isless(a.b, b.b)
end

# ---- Display ----
Base.show(io::IO, p::MetricTensor) = print(io, "g($(p.a),$(p.b))")
Base.show(io::IO, p::Pair{LorentzIndex, Momentum}) = print(io, "$(p.b)^$(p.a)")
Base.show(io::IO, p::Pair{Momentum, Momentum}) = print(io, "$(p.a).$(p.b)")
Base.show(io::IO, p::Pair) = print(io, "Pair($(p.a),$(p.b))")

# ---- Convenience constructors (user-facing API) ----
SP(p::Symbol, q::Symbol) = pair(Momentum(p), Momentum(q))
FV(p::Symbol, mu::Symbol) = pair(LorentzIndex(mu), Momentum(p))
MT(mu::Symbol, nu::Symbol) = pair(LorentzIndex(mu), LorentzIndex(nu))

# D-dimensional variants
SPD(p::Symbol, q::Symbol) = pair(Momentum(p, DimD()), Momentum(q, DimD()))
FVD(p::Symbol, mu::Symbol) = pair(LorentzIndex(mu, DimD()), Momentum(p, DimD()))
MTD(mu::Symbol, nu::Symbol) = pair(LorentzIndex(mu, DimD()), LorentzIndex(nu, DimD()))
