# Feynfeld.jl — ExpandScalarProduct
#
# Bilinear expansion of Pair objects containing MomentumSum.
# Pair(p+q, k) → Pair(p, k) + Pair(q, k)
# Pair(p+q, k+l) → Pair(p,k) + Pair(p,l) + Pair(q,k) + Pair(q,l)
#
# Ref: FeynCalc Lorentz/ExpandScalarProduct.m

export expand_scalar_product

"""
    expand_scalar_product(p::Pair; ctx=SPContext()) -> Vector{Tuple{Rational{Int}, Any}}

Expand a Pair bilinearly over MomentumSum arguments.

Returns a list of `(coefficient, term)` pairs where each term is a `Pair`
(or a scalar from SPContext lookup). Returns `[(1//1, p)]` unchanged if
no MomentumSum is present.

# Examples
```julia
pq = Momentum(:p) + Momentum(:q)
fv = Pair(LorentzIndex(:μ), pq)
expand_scalar_product(fv)
# → [(1//1, FV(:p,:μ)), (1//1, FV(:q,:μ))]

sp = Pair(pq, Momentum(:k))
expand_scalar_product(sp)
# → [(1//1, SP(:k,:p)), (1//1, SP(:k,:q))]
```
"""
function expand_scalar_product(p::Feynfeld.Pair; ctx::SPContext=SPContext())
    a_terms = _expand_slot(p.a)
    b_terms = _expand_slot(p.b)

    # If neither slot was a MomentumSum, return unchanged
    (length(a_terms) == 1 && length(b_terms) == 1 &&
     a_terms[1][1] == 1 // 1 && b_terms[1][1] == 1 // 1) && return [(1 // 1, p)]

    result = Tuple{Rational{Int},Any}[]
    for (ca, aa) in a_terms
        for (cb, bb) in b_terms
            coeff = ca * cb
            r = pair(aa, bb)
            r == 0 && continue
            # Apply SP registry
            if r isa Feynfeld.Pair && r.a isa Momentum && r.b isa Momentum
                val = get_sp(ctx, r.a.name, r.b.name)
                val !== nothing && (r = val)
            end
            push!(result, (coeff, r))
        end
    end

    isempty(result) && return [(0 // 1, 0)]
    result
end

"""Decompose a PairArg slot into (coefficient, atomic PairArg) terms."""
_expand_slot(m::Momentum) = [(1 // 1, m)]
_expand_slot(li::LorentzIndex) = [(1 // 1, li)]
_expand_slot(ms::MomentumSum) = [(c, Momentum(s, ms.dim)) for (c, s) in ms.terms]
