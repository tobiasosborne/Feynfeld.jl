# Feynfeld.jl — Dirac algebra type system
#
# DiracGamma is a tagged struct with three slot variants:
#   LISlot   — gamma^mu (Lorentz index)
#   MomSlot  — p-slash (slashed momentum)
#   SpecialSlot — gamma^5, gamma^6 (P_R), gamma^7 (P_L)
#
# Ported from FeynCalc SharedObjects.m DiracGamma[...].
# BMHV projection applies at construction via dim_contract.

export DiracGamma, dirac_gamma, Spinor, DiracChain
export LISlot, MomSlot, SpecialSlot, DiracSlot, gamma_dim
export GA, GAD, GAE, GS, GSD, GSE, GA5, GA6, GA7

# ── Slot types ───────────────────────────────────────────────────────

"""gamma^mu: carries a LorentzIndex."""
struct LISlot
    index::LorentzIndex
end

"""p-slash = gamma^mu p_mu: carries a Momentum (or MomentumSum)."""
struct MomSlot
    mom::Union{Momentum,MomentumSum}
end

"""gamma^5, gamma^6=(1+g5)/2, gamma^7=(1-g5)/2."""
struct SpecialSlot
    id::Int  # 5, 6, or 7
end

const DiracSlot = Union{LISlot, MomSlot, SpecialSlot}

# ── DiracGamma ───────────────────────────────────────────────────────

"""
    DiracGamma(slot::DiracSlot) <: FeynExpr

A Dirac gamma matrix. The slot determines the variant:
- `LISlot(LorentzIndex)` → gamma^mu
- `MomSlot(Momentum)` → slashed momentum p-slash
- `SpecialSlot(5/6/7)` → gamma^5, P_R, P_L
"""
struct DiracGamma <: FeynExpr
    slot::DiracSlot
end

"""Get the dimension slot of a DiracGamma."""
gamma_dim(g::DiracGamma) = _slot_dim(g.slot)
_slot_dim(s::LISlot) = s.index.dim
_slot_dim(s::MomSlot) = get_dim(s.mom)
_slot_dim(::SpecialSlot) = Dim4()

Base.hash(g::DiracGamma, h::UInt) = hash(g.slot, hash(:DiracGamma, h))
Base.hash(s::LISlot, h::UInt) = hash(s.index, hash(:LISlot, h))
Base.hash(s::MomSlot, h::UInt) = hash(s.mom, hash(:MomSlot, h))
Base.hash(s::SpecialSlot, h::UInt) = hash(s.id, hash(:SpecialSlot, h))

# ── Factory with BMHV projection ────────────────────────────────────

"""
    dirac_gamma(index_or_momentum_or_int, [dim]) -> Union{DiracGamma, Int}

Create a DiracGamma with BMHV projection. Returns 0 when BMHV vanishes.
"""
function dirac_gamma(li::LorentzIndex)
    DiracGamma(LISlot(li))
end

function dirac_gamma(li::LorentzIndex, dim::DimSlot)
    if !(li.dim === dim)
        dp::Union{DimSlot,Nothing} = dim_contract(li.dim, dim)
        dp === nothing && return 0
        li = LorentzIndex(li.name, dp::DimSlot)
    end
    DiracGamma(LISlot(li))
end

function dirac_gamma(m::Momentum)
    DiracGamma(MomSlot(m))
end

function dirac_gamma(m::Momentum, dim::DimSlot)
    if !(m.dim === dim)
        dp::Union{DimSlot,Nothing} = dim_contract(m.dim, dim)
        dp === nothing && return 0
        m = Momentum(m.name, dp::DimSlot)
    end
    DiracGamma(MomSlot(m))
end

function dirac_gamma(ms::MomentumSum)
    DiracGamma(MomSlot(ms))
end

function dirac_gamma(id::Int)
    id in (5, 6, 7) || error("DiracGamma special id must be 5, 6, or 7; got $id")
    DiracGamma(SpecialSlot(id))
end

# ── Convenience constructors ─────────────────────────────────────────

"""GA(μ) — 4D gamma^μ."""
GA(mu::Symbol) = dirac_gamma(LorentzIndex(mu))
"""GAD(μ) — D-dimensional gamma^μ."""
GAD(mu::Symbol) = dirac_gamma(LorentzIndex(mu, DimD()))
"""GAE(μ) — evanescent (D-4)-dim gamma^μ."""
GAE(mu::Symbol) = dirac_gamma(LorentzIndex(mu, DimDm4()))
"""GS(p) — 4D slashed momentum p-slash."""
GS(p::Symbol) = dirac_gamma(Momentum(p))
"""GSD(p) — D-dimensional slashed momentum."""
GSD(p::Symbol) = dirac_gamma(Momentum(p, DimD()))
"""GSE(p) — evanescent slashed momentum."""
GSE(p::Symbol) = dirac_gamma(Momentum(p, DimDm4()))
"""GA5() — gamma^5."""
GA5() = dirac_gamma(5)
"""GA6() — chiral projector (1+gamma^5)/2."""
GA6() = dirac_gamma(6)
"""GA7() — chiral projector (1-gamma^5)/2."""
GA7() = dirac_gamma(7)

# ── Spinor ───────────────────────────────────────────────────────────

"""
    Spinor(kind, momentum, mass)

A Dirac spinor. `kind` is `:u`, `:v`, `:ubar`, or `:vbar`.

Ref: FeynCalc SharedObjects.m Spinor[Momentum[p], m, 1].
"""
struct Spinor <: FeynExpr
    kind::Symbol
    momentum::Momentum
    mass::Any
    function Spinor(kind::Symbol, mom::Momentum, mass=0)
        kind in (:u, :v, :ubar, :vbar) || error("Spinor kind must be :u, :v, :ubar, or :vbar")
        new(kind, mom, mass)
    end
end

Base.hash(s::Spinor, h::UInt) = hash(s.mass, hash(s.momentum, hash(s.kind, hash(:Spinor, h))))

# ── DiracChain (minimal placeholder for Phase 1b.2) ─────────────────

"""
    DiracChain(elements)

A non-commutative chain of DiracGamma and Spinor objects.
Full DOT semantics (scalar extraction, expansion) in Phase 1b.2.
"""
struct DiracChain <: FeynExpr
    elements::Vector{Union{DiracGamma,Spinor}}
end

Base.length(c::DiracChain) = length(c.elements)
