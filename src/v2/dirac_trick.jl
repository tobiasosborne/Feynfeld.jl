# DiracTrick: simplify gamma chains by contracting repeated Lorentz indices.
# This is NOT trace — it operates on matrix-valued DiracExpr.
#
# Key identities in D dimensions:
#   γ^μ γ_μ = D                                       (n=0)
#   γ^μ γ^a γ_μ = -(D-2) γ^a                          (n=1)
#   γ^μ γ^a γ^b γ_μ = 4 g^{ab} + (D-4) γ^a γ^b       (n=2)
#
# For the 1-loop self-energy, n=0 and n=1 suffice.

# Contract repeated Lorentz indices in a DiracExpr → DiracExpr
function dirac_trick(de::DiracExpr)
    result = DiracExpr()
    for (coeff, chain) in de.terms
        simplified = _trick_chain(coeff, chain)
        result = result + simplified
    end
    simplify(result)
end

# Simplify a single chain: find and contract repeated Lorentz indices
function _trick_chain(coeff::AlgSum, chain::DiracChain)
    gs = DiracGamma[g for g in chain.elements if g isa DiracGamma]
    _trick_gammas(coeff, gs)
end

function _trick_gammas(coeff::AlgSum, gs::Vector{DiracGamma})
    # First expand any MomSumSlot
    idx = findfirst(g -> g isa DiracGamma{MomSumSlot}, gs)
    if idx !== nothing
        g = gs[idx]::DiracGamma{MomSumSlot}
        result = DiracExpr()
        for (c, m) in g.slot.mom.terms
            expanded = copy(gs)
            expanded[idx] = GS(m)
            result = result + _trick_gammas(c * coeff, expanded)
        end
        return result
    end

    # Find repeated Lorentz index
    for (i, gi) in enumerate(gs)
        li = lorentz_index(gi)
        li === nothing && continue
        for j in (i+1):length(gs)
            lj = lorentz_index(gs[j])
            lj === nothing && continue
            li == lj || continue

            # Found γ^μ ... γ_μ at positions i and j
            n = j - i - 1  # number of gammas between
            inner = gs[i+1:j-1]
            left = gs[1:i-1]
            right = gs[j+1:end]
            dim = li.dim

            return _apply_sandwich(coeff, left, inner, right, dim)
        end
    end

    # No repeated index found — return as-is
    DiracExpr([(coeff, DiracChain(gs))])
end

# Apply the sandwich identity: γ^μ [inner] γ_μ
function _apply_sandwich(coeff, left, inner, right, dim)
    n = length(inner)
    d = dim_trace(dim)  # D or 4, as Coeff

    if n == 0
        # γ^μ γ_μ = D → scalar D times identity
        # Result: D * coeff * [left . right]
        new_gs = vcat(left, right)
        new_coeff = isempty(new_gs) ? d * coeff : d * coeff
        return _trick_gammas(d * coeff, DiracGamma[new_gs...])
    end

    if n == 1
        # γ^μ γ^a γ_μ = -(D-2) γ^a
        # Result: (2-D) * coeff * [left . γ^a . right]
        factor = add_coeff(2//1, isa(d, DimPoly) ? -d : -d)
        factor = normalise_coeff(factor)
        new_gs = DiracGamma[left; inner; right]
        return _trick_gammas(factor * coeff, new_gs)
    end

    if n == 2
        # γ^μ γ^a γ^b γ_μ = 4 g^{ab} + (D-4) γ^a γ^b
        # Term 1: 4 * g^{ab} * [left . right]  (scalar, no gammas from inner)
        # Term 2: (D-4) * [left . γ^a . γ^b . right]
        ga, gb = inner[1], inner[2]
        p = gamma_pair(ga, gb)

        # Term 1: metric contraction
        term1 = if p !== nothing && !(p isa Number && iszero(p))
            p_coeff = p isa Number ? p * coeff : alg(p) * coeff
            _trick_gammas(4 * p_coeff, DiracGamma[left; right])
        else
            DiracExpr()
        end

        # Term 2: (D-4) times original inner
        dm4 = add_coeff(d, -4//1)
        dm4 = normalise_coeff(dm4)
        term2 = if !_coeff_iszero(dm4)
            new_gs = DiracGamma[left; inner; right]
            _trick_gammas(dm4 * coeff, new_gs)
        else
            DiracExpr()  # vanishes in 4 dimensions
        end

        return term1 + term2
    end

    # n ≥ 3: not implemented yet (would need general Mertig/Boehm/Denner formula)
    error("DiracTrick: γ^μ [$(n) gammas] γ_μ not implemented for n ≥ 3")
end
