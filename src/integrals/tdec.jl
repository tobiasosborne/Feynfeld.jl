# Feynfeld.jl — Tensor decomposition (Tdec)
#
# Decomposes tensor loop integrals into scalar integrals times
# tensor structures built from external momenta and g^{μν}.
#
# For a rank-R integral with loop momentum q and external momenta p_i:
#   ∫ q^{μ1}...q^{μR} / denom = Σ C_i * (tensor structure)_i
#
# Phase 1d: rank 1 and 2 implemented. Higher ranks deferred.
#
# Ref: FeynCalc LoopIntegrals/Tdec.m, Passarino-Veltman (1979)

export tdec

"""
    tdec(indices, external_momenta; dim=DimD()) -> Vector{Tuple{Any, Any}}

Tensor decomposition of a loop integral.

`indices`: list of `(loop_momentum, lorentz_index)` pairs specifying the tensor rank.
`external_momenta`: list of external momentum symbols.

Returns a list of `(PaVe_coefficient, tensor_structure)` pairs.

# Examples
```julia
# Rank 1: ∫ q^μ / (q²-m²)((q-p)²-m²) = B1 * p^μ
tdec([(:q, :μ)], [:p])
# → [(:B1, FVD(:p, :μ))]
```
"""
function tdec(indices::Vector, ext_mom::Vector{Symbol}; dim::DimSlot=DimD())
    rank = length(indices)
    n_ext = length(ext_mom)

    rank == 0 && return [(1, 1)]  # scalar: trivially 1 * scalar_integral

    if rank == 1
        return _tdec_rank1(indices[1], ext_mom, dim)
    elseif rank == 2
        return _tdec_rank2(indices, ext_mom, dim)
    else
        error("Tdec for rank $rank not yet implemented (Phase 1d supports rank ≤ 2)")
    end
end

"""Rank-1 decomposition: ∫ q^μ = Σ_i C_i p_i^μ."""
function _tdec_rank1(idx_pair, ext_mom, dim)
    _, mu = idx_pair
    mu_li = LorentzIndex(mu, dim)
    # Each external momentum contributes one term: C_i * p_i^μ
    result = Tuple{Symbol,Any}[]
    for (i, p) in enumerate(ext_mom)
        coeff = Symbol("C", i)  # placeholder coefficient label
        push!(result, (coeff, Feynfeld.Pair(mu_li, Momentum(p, dim))))
    end
    result
end

"""Rank-2 decomposition: ∫ q^μ q^ν = C_00 g^{μν} + Σ_{ij} C_ij p_i^μ p_j^ν."""
function _tdec_rank2(indices, ext_mom, dim)
    _, mu = indices[1]
    _, nu = indices[2]
    mu_li = LorentzIndex(mu, dim)
    nu_li = LorentzIndex(nu, dim)

    result = Tuple{Symbol,Any}[]

    # Metric term: C_00 * g^{μν}
    push!(result, (:C00, Feynfeld.Pair(mu_li, nu_li)))

    # External momentum terms: C_ij * p_i^μ * p_j^ν
    for (i, pi) in enumerate(ext_mom)
        for (j, pj) in enumerate(ext_mom)
            j < i && continue  # C_ij = C_ji by symmetry
            coeff = Symbol("C", i, j)
            fv_mu = Feynfeld.Pair(mu_li, Momentum(pi, dim))
            fv_nu = Feynfeld.Pair(nu_li, Momentum(pj, dim))
            push!(result, (coeff, (fv_mu, fv_nu)))
        end
    end
    result
end
