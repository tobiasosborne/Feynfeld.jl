# Feynfeld.jl — Pair: the universal Lorentz bilinear
#
# Pair is the single most important type in the Lorentz algebra.
# Every Lorentz-covariant scalar, vector, and tensor reduces to Pair.
#
# Ported from FeynCalc SharedObjects.m Pair[x, y].
# Key properties:
#   - Symmetric (Orderless): Pair(a,b) == Pair(b,a)
#   - BMHV projection at construction: mismatched dims are projected
#   - Meaning depends on argument types:
#       Pair(LorentzIndex, LorentzIndex)  = metric tensor g^{μν}
#       Pair(Momentum, Momentum)          = scalar product p·q
#       Pair(LorentzIndex, Momentum)      = four-vector component p^μ

# Note: Pair is NOT exported (conflicts with Base.Pair).
# Users access it as Feynfeld.Pair or via convenience constructors.
export pair
export SP, SPD, SPE, FV, FVD, FVE, MT, MTD, MTE

# ── Helpers ──────────────────────────────────────────────────────────

"""Get the DimSlot of a PairArg."""
get_dim(li::LorentzIndex) = li.dim
get_dim(m::Momentum) = m.dim
get_dim(ms::MomentumSum) = ms.dim

"""Return a copy with a new DimSlot."""
set_dim(li::LorentzIndex, d::DimSlot) = LorentzIndex(li.name, d)
set_dim(m::Momentum, d::DimSlot) = Momentum(m.name, d)
set_dim(ms::MomentumSum, d::DimSlot) = MomentumSum(ms.terms, d)

"""Check if a PairArg represents zero momentum."""
_is_zero(::LorentzIndex) = false
_is_zero(::Momentum) = false
_is_zero(ms::MomentumSum) = isempty(ms.terms)

# Type ordering for canonical form: LorentzIndex < Momentum < MomentumSum
_pair_type_order(::LorentzIndex) = 1
_pair_type_order(::Momentum) = 2
_pair_type_order(::MomentumSum) = 3

# ── Pair struct ──────────────────────────────────────────────────────

"""
    Pair(a::PairArg, b::PairArg) <: FeynExpr

The universal Lorentz bilinear. Symmetric: `Pair(a,b) == Pair(b,a)`.

Construction enforces:
1. BMHV dimension projection (mismatched dims → projected via dim_contract)
2. Canonical ordering (LorentzIndex before Momentum, then alphabetical)

Errors if BMHV projection vanishes (4 ∩ (D-4) = 0).
Use `pair()` when vanishing is expected (returns 0 instead of erroring).
"""
struct Pair <: FeynExpr
    a::PairArg
    b::PairArg
    function Pair(x::PairArg, y::PairArg)
        # BMHV projection
        dx, dy = get_dim(x), get_dim(y)
        if !(dx === dy)
            dp::Union{DimSlot,Nothing} = dim_contract(dx, dy)
            dp === nothing && error(
                "Pair vanishes: $(typeof(dx)) ∩ $(typeof(dy)) = 0. " *
                "Use pair() for expressions that may vanish."
            )
            x = set_dim(x, dp::DimSlot)
            y = set_dim(y, dp::DimSlot)
        end
        # Canonical ordering (Orderless)
        ox, oy = _pair_type_order(x), _pair_type_order(y)
        if ox > oy || (ox == oy && isless(y, x))
            x, y = y, x
        end
        new(x, y)
    end
end

Base.hash(p::Pair, h::UInt) = hash(p.b, hash(p.a, hash(:Pair, h)))

# ── Factory function ─────────────────────────────────────────────────

"""
    pair(a::PairArg, b::PairArg) -> Union{Pair, Int}

Create a Pair with BMHV projection. Returns `0` when the projection
vanishes (4 ∩ (D-4) = 0), instead of erroring.

Ref: FeynCalc SharedObjects.m lines 2361–2381 (Pair DownValues).
"""
function pair(x::PairArg, y::PairArg)
    # Zero momentum check (empty MomentumSum)
    (_is_zero(x) || _is_zero(y)) && return 0
    # BMHV projection
    dx, dy = get_dim(x), get_dim(y)
    if !(dx === dy)
        dp::Union{DimSlot,Nothing} = dim_contract(dx, dy)
        dp === nothing && return 0
        x = set_dim(x, dp::DimSlot)
        y = set_dim(y, dp::DimSlot)
    end
    Pair(x, y)
end

# ── Convenience constructors ─────────────────────────────────────────
# These mirror FeynCalc's FCE shorthand: SP, FV, MT etc.

# Metric tensor g^{μν}
"""MT(μ, ν) — 4-dimensional metric tensor g^{μν}."""
MT(mu::Symbol, nu::Symbol) = Pair(LorentzIndex(mu), LorentzIndex(nu))
"""MTD(μ, ν) — D-dimensional metric tensor g^{μν}."""
MTD(mu::Symbol, nu::Symbol) = Pair(LorentzIndex(mu, DimD()), LorentzIndex(nu, DimD()))
"""MTE(μ, ν) — evanescent (D-4)-dimensional metric tensor."""
MTE(mu::Symbol, nu::Symbol) = Pair(LorentzIndex(mu, DimDm4()), LorentzIndex(nu, DimDm4()))

# Four-vector component p^μ
"""FV(p, μ) — 4-dimensional four-vector component p^μ."""
FV(p::Symbol, mu::Symbol) = Pair(LorentzIndex(mu), Momentum(p))
"""FVD(p, μ) — D-dimensional four-vector component p^μ."""
FVD(p::Symbol, mu::Symbol) = Pair(LorentzIndex(mu, DimD()), Momentum(p, DimD()))
"""FVE(p, μ) — evanescent four-vector component."""
FVE(p::Symbol, mu::Symbol) = Pair(LorentzIndex(mu, DimDm4()), Momentum(p, DimDm4()))

# Scalar product p·q
"""SP(p, q) — 4-dimensional scalar product p·q."""
SP(p::Symbol, q::Symbol) = Pair(Momentum(p), Momentum(q))
"""SP(p) — 4-dimensional p²."""
SP(p::Symbol) = SP(p, p)
"""SPD(p, q) — D-dimensional scalar product p·q."""
SPD(p::Symbol, q::Symbol) = Pair(Momentum(p, DimD()), Momentum(q, DimD()))
"""SPD(p) — D-dimensional p²."""
SPD(p::Symbol) = SPD(p, p)
"""SPE(p, q) — evanescent scalar product."""
SPE(p::Symbol, q::Symbol) = Pair(Momentum(p, DimDm4()), Momentum(q, DimDm4()))
"""SPE(p) — evanescent p²."""
SPE(p::Symbol) = SPE(p, p)
