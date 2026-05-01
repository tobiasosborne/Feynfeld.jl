# Levi-Civita tensor contraction: ε·ε → determinant of metric tensors.
#
# Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.21)
# "ε^{μνρσ} ε_{μ'ν'ρ'σ'} = -det|g^μ_{μ'} g^μ_{ν'} g^μ_{ρ'} g^μ_{σ'}|"
#                                 |g^ν_{μ'} g^ν_{ν'} ...              |
#                                 |...                                |
# Cross-check: refs/FeynCalc/Tests/Lorentz/EpsContract.test (41 tests)
#
# Convention: sign² = 1 (same as FeynCalc $LeviCivitaSign = -1, (-1)² = 1).
# Overall: ε₁ · ε₂ = -det[pair(a_i, b_j)]_{4×4}

# ---- S₄ permutations with parity (precomputed) ----
const _S4_PERMS = let
    perms = Tuple{NTuple{4,Int}, Int}[]
    for p in [[1,2,3,4], [1,2,4,3], [1,3,2,4], [1,3,4,2], [1,4,2,3], [1,4,3,2],
              [2,1,3,4], [2,1,4,3], [2,3,1,4], [2,3,4,1], [2,4,1,3], [2,4,3,1],
              [3,1,2,4], [3,1,4,2], [3,2,1,4], [3,2,4,1], [3,4,1,2], [3,4,2,1],
              [4,1,2,3], [4,1,3,2], [4,2,1,3], [4,2,3,1], [4,3,1,2], [4,3,2,1]]
        # Count inversions for parity
        inv = sum(p[i] > p[j] for i in 1:4 for j in (i+1):4)
        sign = iseven(inv) ? 1 : -1
        push!(perms, (Tuple(p), sign))
    end
    perms
end

"""
    eps_contract(s::AlgSum) → AlgSum

Pre-pass: find pairs of Eps factors in each term, replace with the
determinant formula. Call BEFORE the per-index contraction loop.
"""
function eps_contract(s::AlgSum)
    result = AlgSum()
    for (fk, c) in s.terms
        contracted = _eps_contract_term(fk.factors, c)
        result = result + contracted
    end
    result
end

function _eps_contract_term(factors::Vector{AlgFactor}, coeff)
    # Find Eps factors
    eps_idx = [i for (i, f) in enumerate(factors) if f isa Eps]
    length(eps_idx) < 2 && return alg_from_factors(factors, coeff)

    # Contract first pair found (contract() will re-enter for additional pairs)
    i, j = eps_idx[1], eps_idx[2]
    e1, e2 = factors[i]::Eps, factors[j]::Eps

    # Remaining factors (everything except the two Eps)
    remaining = AlgFactor[factors[k] for k in eachindex(factors) if k != i && k != j]

    # Compute -det[pair(a_i, b_j)]_{4×4}
    args1 = PairArg[e1.a, e1.b, e1.c, e1.d]
    args2 = PairArg[e2.a, e2.b, e2.c, e2.d]
    det_sum = _det4x4_pairs(args1, args2)

    # Multiply det by -1 (the overall sign) and by remaining factors + coefficient
    remaining_alg = alg_from_factors(remaining, coeff)
    -det_sum * remaining_alg
end

"""
    _det4x4_pairs(a, b) → AlgSum

Symbolic 4×4 determinant where entry (i,j) = pair(a[i], b[j]).
Uses Leibniz formula: det = Σ_σ sgn(σ) Π_i M_{i,σ(i)}
"""
function _det4x4_pairs(a::Vector{PairArg}, b::Vector{PairArg})
    result = AlgSum()
    for (perm, sign) in _S4_PERMS
        term = alg(sign)
        for i in 1:4
            ai, bi = a[i], b[perm[i]]
            # BMHV vanishing for mixed-dim LI×LI: skip this permutation term
            if ai isa LorentzIndex && bi isa LorentzIndex &&
               dim_contract(ai.dim, bi.dim) === nothing
                term = AlgSum()
                break
            end
            term = term * alg(pair(ai, bi))
        end
        result = result + term
    end
    result
end
