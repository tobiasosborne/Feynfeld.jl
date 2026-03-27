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
            other_factors = [factors[j] for j in eachindex(factors) if j != i]
            result = AlgSum()
            for (c, new_f) in expanded
                new_factors = isnothing(new_f) ? copy(other_factors) : [other_factors; new_f]
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
_try_expand(::Eps) = nothing

# MomentumSum in first slot
function _try_expand(p::Pair{MomentumSum, A}) where {A<:PairArg}
    result = Tuple{Any, Union{AlgFactor, Nothing}}[]
    for (c, m) in p.a.terms
        new_p = pair(m, p.b)
        if new_p isa Number
            push!(result, (mul_coeff(c, new_p), nothing))
        else
            push!(result, (c, new_p))
        end
    end
    result
end

# MomentumSum in second slot
function _try_expand(p::Pair{A, MomentumSum}) where {A<:PairArg}
    result = Tuple{Any, Union{AlgFactor, Nothing}}[]
    for (c, m) in p.b.terms
        new_p = pair(p.a, m)
        if new_p isa Number
            push!(result, (mul_coeff(c, new_p), nothing))
        else
            push!(result, (c, new_p))
        end
    end
    result
end

# Both slots are MomentumSum
function _try_expand(p::Pair{MomentumSum, MomentumSum})
    result = Tuple{Any, Union{AlgFactor, Nothing}}[]
    for (c1, m1) in p.a.terms
        for (c2, m2) in p.b.terms
            c = mul_coeff(c1, c2)
            new_p = pair(m1, m2)
            if new_p isa Number
                push!(result, (mul_coeff(c, new_p), nothing))
            else
                push!(result, (c, new_p))
            end
        end
    end
    result
end
