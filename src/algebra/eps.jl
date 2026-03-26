# Feynfeld.jl — Levi-Civita tensor and contraction
#
# Eps represents ε^{μνρσ}, the totally antisymmetric Levi-Civita tensor.
# Args can be LorentzIndex or Momentum (4 slots).
#
# EpsContract: ε^{a1..a4} ε_{b1..b4} = -det[g^{ai,bj}]
# Implements n-shared-index specialisations for 4D and D-dim.
#
# Ref: FeynCalc Lorentz/Eps.m, EpsContract.m, EpsEvaluate.m

export Eps, levi_civita, LC, LCD, eps_contract

"""Default FeynCalc convention: \$LeviCivitaSign = -1."""
const LEVI_CIVITA_SIGN = -1

# ── Eps type ─────────────────────────────────────────────────────────

"""
    Eps(args::NTuple{4, PairArg})

Levi-Civita tensor ε^{a b c d}. Args are LorentzIndex or Momentum.
Construct via `levi_civita()` or `LC()`/`LCD()` convenience functions.
"""
struct Eps <: FeynExpr
    args::NTuple{4,PairArg}
end

"""
    levi_civita(a, b, c, d) -> Union{Eps, Int}

Create an Eps tensor. Returns 0 if:
- Any argument is in the (D-4) evanescent hat space (BMHV)
- Any two arguments are identical (antisymmetry)
"""
function levi_civita(a::PairArg, b::PairArg, c::PairArg, d::PairArg)
    args = (a, b, c, d)
    for x in args
        get_dim(x) isa DimDm4 && return 0
    end
    for i in 1:4, j in (i + 1):4
        args[i] == args[j] && return 0
    end
    Eps(args)
end

Base.hash(e::Eps, h::UInt) = hash(e.args, hash(:Eps, h))

# ── Convenience constructors ─────────────────────────────────────────

"""LC(a, b, c, d) — 4-dimensional Levi-Civita ε^{abcd}."""
function LC(a::Symbol, b::Symbol, c::Symbol, d::Symbol)
    levi_civita(LorentzIndex(a), LorentzIndex(b), LorentzIndex(c), LorentzIndex(d))
end

"""LCD(a, b, c, d) — D-dimensional Levi-Civita ε^{abcd}."""
function LCD(a::Symbol, b::Symbol, c::Symbol, d::Symbol)
    levi_civita(LorentzIndex(a, DimD()), LorentzIndex(b, DimD()),
                LorentzIndex(c, DimD()), LorentzIndex(d, DimD()))
end

# ── EpsContract ──────────────────────────────────────────────────────

"""
    eps_contract(e1::Eps, e2::Eps) -> result

Contract two Levi-Civita tensors using the determinant identity:
  ε^{a1..a4} ε_{b1..b4} = -(\$LeviCivitaSign)² det[g^{ai,bj}]

Currently handles n shared LorentzIndex names (n = 2, 3, 4).
Returns a scalar, Pair, or sum-of-products depending on shared count.

# Examples
```julia
eps_contract(LC(:a,:b,:c,:d), LC(:a,:b,:c,:d))   # → -24
eps_contract(LC(:a,:b,:c,:μ), LC(:a,:b,:c,:ν))   # → -6 * MT(:μ,:ν)
```
"""
function eps_contract(e1::Eps, e2::Eps)
    n, free1, free2 = _eps_shared(e1, e2)
    dim = get_dim(e1.args[1])

    if n == 4
        return _eps_scalar(dim, 4)
    elseif n == 3
        g = pair(free1[1], free2[1])
        g == 0 && return 0
        return (_eps_scalar(dim, 3), g)  # (coeff, metric_pair)
    elseif n == 2
        # -2(g_{a1b1}g_{a2b2} - g_{a1b2}g_{a2b1}) in 4D
        c = _eps_scalar(dim, 2)
        g11 = pair(free1[1], free2[1])
        g22 = pair(free1[2], free2[2])
        g12 = pair(free1[1], free2[2])
        g21 = pair(free1[2], free2[1])
        # Return as list of (coeff, factors) terms
        terms = Tuple{Any,Vector}[]
        _push_product!(terms, c, g11, g22)
        _push_product!(terms, _mul_coeff(-1, c), g12, g21)
        return terms
    else
        error("eps_contract with $n shared indices: general case not yet implemented")
    end
end

# ── Internals ────────────────────────────────────────────────────────

"""Count shared LorentzIndex names between two Eps, return (count, free1, free2)."""
function _eps_shared(e1::Eps, e2::Eps)
    free2 = collect(PairArg, e2.args)
    free1 = PairArg[]
    shared = 0
    for a in e1.args
        matched = false
        if a isa LorentzIndex
            for j in eachindex(free2)
                if free2[j] isa LorentzIndex && a.name == free2[j].name
                    shared += 1
                    deleteat!(free2, j)
                    matched = true
                    break
                end
            end
        end
        matched || push!(free1, a)
    end
    (shared, free1, free2)
end

"""Scalar factor for n-shared eps contraction: -n! in 4D, falling factorial in D-dim."""
function _eps_scalar(::Dim4, n)
    -(LEVI_CIVITA_SIGN)^2 * factorial(n)
end

function _eps_scalar(::DimD, n)
    sign = -(LEVI_CIVITA_SIGN)^2
    if n == 4;     :($sign * D * (D - 1) * (D - 2) * (D - 3))
    elseif n == 3; :($sign * (D - 1) * (D - 2) * (D - 3))
    elseif n == 2; :($sign * (D - 2) * (D - 3))
    else error("D-dim eps scalar for n=$n not implemented")
    end
end

"""Push a product term (coeff, [g1, g2]) into terms, handling zeros."""
function _push_product!(terms, coeff, g1, g2)
    (coeff == 0 || g1 == 0 || g2 == 0) && return
    factors = Any[]
    g1 isa Feynfeld.Pair && push!(factors, g1)
    g2 isa Feynfeld.Pair && push!(factors, g2)
    push!(terms, (coeff, factors))
end
