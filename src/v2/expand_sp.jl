# Expand scalar products of MomentumSums bilinearly.
# Dispatch on parametric Pair types — clean and type-stable.

function expand_scalar_product(s::AlgSum)
    result = AlgSum()
    for (fk, c) in s.terms
        expanded = _expand_term(fk.factors, c)
        result = result + expanded
    end
    result
end

function _expand_term(factors::Vector{AlgFactor}, coeff)
    # Find first expandable factor
    for (i, f) in enumerate(factors)
        expanded = _try_expand(f)
        if expanded !== nothing
            # Replace factor i with expanded terms, recurse
            other_factors = AlgFactor[factors[j] for j in eachindex(factors) if j != i]
            result = AlgSum()
            for (c, new_f) in expanded
                new_factors = AlgFactor[other_factors; new_f]
                sub = _expand_term(new_factors, mul_coeff(coeff, c))
                result = result + sub
            end
            return result
        end
    end
    # Nothing to expand
    alg_from_factors(factors, coeff)
end

# Dispatch on Pair type: only expand if a MomentumSum is involved
_try_expand(::Pair{Momentum, Momentum}) = nothing
_try_expand(::Pair{LorentzIndex, Momentum}) = nothing
_try_expand(::Pair{LorentzIndex, LorentzIndex}) = nothing
# Eps with MomentumSum slot: expand by linearity ε(αK₁+βK₂, b, c, d) = α ε(K₁,...) + β ε(K₂,...)
function _try_expand(e::Eps)
    for (slot_idx, slot) in enumerate((e.a, e.b, e.c, e.d))
        slot isa MomentumSum || continue
        result = Tuple{Coeff, AlgFactor}[]
        for (c, m) in slot.terms
            slots = PairArg[e.a, e.b, e.c, e.d]
            slots[slot_idx] = m
            push!(result, (c, Eps(slots[1], slots[2], slots[3], slots[4])))
        end
        return result
    end
    nothing  # no MomentumSum slots
end

# MomentumSum in first slot
function _try_expand(p::Pair{MomentumSum, A}) where {A<:PairArg}
    result = Tuple{Coeff, AlgFactor}[]
    for (c, m) in p.a.terms
        push!(result, (c, pair(m, p.b)))
    end
    result
end

# MomentumSum in second slot
function _try_expand(p::Pair{A, MomentumSum}) where {A<:PairArg}
    result = Tuple{Coeff, AlgFactor}[]
    for (c, m) in p.b.terms
        push!(result, (c, pair(p.a, m)))
    end
    result
end

# Both slots are MomentumSum
function _try_expand(p::Pair{MomentumSum, MomentumSum})
    result = Tuple{Coeff, AlgFactor}[]
    for (c1, m1) in p.a.terms
        for (c2, m2) in p.b.terms
            push!(result, (mul_coeff(c1, c2), pair(m1, m2)))
        end
    end
    result
end
