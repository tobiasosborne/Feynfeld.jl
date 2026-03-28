# DiracTrick: simplify gamma chains by contracting repeated Lorentz indices.
# This is NOT trace — it operates on matrix-valued DiracExpr.
#
# Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9)
# "γ^μ Γ^(l) γ_μ = (-1)^l {(D-2l)Γ^(l) - 4 Σ_{i<j} (-1)^{j-i} Γ_{ij}^(l) g_{μ_i μ_j}}"
#
# Explicit cases (n = number of gammas between contracted pair):
#   n=0: γ^μ γ_μ = D                                        [Eq. (2.2)]
#   n=1: γ^μ γ^a γ_μ = (2-D) γ^a                           [Eq. (2.9), l=1]
#   n=2: γ^μ γ^a γ^b γ_μ = 4g^{ab} + (D-4) γ^a γ^b        [Eq. (2.9), l=2]
#   n=3: γ^μ γ^a γ^b γ^c γ_μ = -(D-4) γ^a γ^b γ^c - 2 γ^c γ^b γ^a
#   n=4: γ^μ γ^a γ^b γ^c γ^d γ_μ = (D-4) γ^a γ^b γ^c γ^d
#         + 2 γ^c γ^b γ^a γ^d + 2 γ^d γ^a γ^b γ^c
#   n≥5: general Eq. (2.9)

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

    if n == 3
        # Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9), l=3
        # "γ^μ γ^a γ^b γ^c γ_μ = -(D-4) γ^a γ^b γ^c - 2 γ^c γ^b γ^a"
        # Verified: refs/FeynCalc/Tests/Dirac/DiracTrick.test, fcstDiracTrickThreeFreeIndices-ID1
        a, b, c = inner[1], inner[2], inner[3]
        dm4 = add_coeff(d, -4//1) |> normalise_coeff

        # Term 1: -(D-4) γ^a γ^b γ^c
        term1 = if !_coeff_iszero(dm4)
            _trick_gammas(-dm4 * coeff, DiracGamma[left; a; b; c; right])
        else
            DiracExpr()
        end

        # Term 2: -2 γ^c γ^b γ^a
        term2 = _trick_gammas(-2 * coeff, DiracGamma[left; c; b; a; right])

        return term1 + term2
    end

    if n == 4
        # Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9), l=4
        # "γ^μ γ^a γ^b γ^c γ^d γ_μ = (D-4) γ^a γ^b γ^c γ^d
        #   + 2 γ^c γ^b γ^a γ^d + 2 γ^d γ^a γ^b γ^c"
        # Verified: refs/FeynCalc/Tests/Dirac/DiracTrick.test, fcstDiracTrickFourFreeIndices-ID1
        a, b, c, e = inner[1], inner[2], inner[3], inner[4]
        dm4 = add_coeff(d, -4//1) |> normalise_coeff

        term1 = if !_coeff_iszero(dm4)
            _trick_gammas(dm4 * coeff, DiracGamma[left; a; b; c; e; right])
        else
            DiracExpr()
        end
        term2 = _trick_gammas(2 * coeff, DiracGamma[left; c; b; a; e; right])
        term3 = _trick_gammas(2 * coeff, DiracGamma[left; e; a; b; c; right])

        return term1 + term2 + term3
    end

    # n ≥ 5: general formula
    # Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9)
    # "γ^μ Γ^(l) γ_μ = (-1)^l {(D-2l)Γ^(l) - 4 Σ_{i<j} (-1)^{j-i} Γ_{ij}^(l) g_{μ_i μ_j}}"
    sign_l = iseven(n) ? 1 : -1
    d_minus_2l = add_coeff(d, Rational{Int}(-2 * n)) |> normalise_coeff

    # First term: (-1)^l (D-2l) Γ^(l) [left . inner . right]
    result = if !_coeff_iszero(d_minus_2l)
        _trick_gammas(sign_l * d_minus_2l * coeff, DiracGamma[left; inner; right])
    else
        DiracExpr()
    end

    # Second term: -4(-1)^l Σ_{i<j} (-1)^{j-i} g_{μ_i μ_j} Γ_{ij}^(l)
    overall = -4 * sign_l
    for i in 1:n-1
        for j in i+1:n
            pair_sign = iseven(j - i) ? 1 : -1
            p = gamma_pair(inner[i], inner[j])
            p === nothing && continue

            remaining = DiracGamma[inner[k] for k in 1:n if k != i && k != j]
            chain = DiracGamma[left; remaining; right]
            p_coeff = p isa Number ? (overall * pair_sign * p) * coeff :
                                    (overall * pair_sign) * (alg(p) * coeff)
            result = result + _trick_gammas(p_coeff, chain)
        end
    end

    return result
end
