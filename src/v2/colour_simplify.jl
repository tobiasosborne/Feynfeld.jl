# Colour simplification: delta contraction, structure constant identities.
# All operations return AlgSum or modify AlgSum in place.

# ---- Contract colour deltas in an AlgSum ----
# δ^{ab} contracts with any factor sharing index a or b, replacing it.
function contract_colour(s::AlgSum; N::Int=3)
    result = AlgSum()
    for (fk, c) in s.terms
        contracted = _contract_colour_term(fk.factors, c, N)
        result = result + contracted
    end
    result
end

function _contract_colour_term(factors::Vector{AlgFactor}, coeff, N::Int)
    factors = copy(factors)
    coeff_acc = coeff
    changed = true
    while changed
        changed = false

        # Look for delta self-contractions (δ^{aa})
        for (i, f) in enumerate(factors)
            if f isa SUNDelta && f.a == f.b
                coeff_acc = mul_coeff(coeff_acc, Rational{Int}(N^2 - 1))
                deleteat!(factors, i)
                changed = true
                break
            elseif f isa FundDelta && f.i == f.j
                coeff_acc = mul_coeff(coeff_acc, Rational{Int}(N))
                deleteat!(factors, i)
                changed = true
                break
            end
        end
        changed && continue

        # Look for δ^{ab} that shares an index with another factor
        for (i, f) in enumerate(factors)
            f isa SUNDelta || continue
            for (j, g) in enumerate(factors)
                i == j && continue  # same position, not same object
                replaced = _delta_substitute(f, g)
                if replaced !== nothing
                    hi, lo = max(i,j), min(i,j)
                    deleteat!(factors, hi)
                    deleteat!(factors, lo)
                    if replaced !== :removed
                        push!(factors, replaced)
                    end
                    changed = true
                    break
                end
            end
            changed && break
        end
        changed && continue

        # f^{acd} f^{bcd} = N δ^{ab} (two shared adjoint indices)
        for i in 1:length(factors)
            factors[i] isa SUNF || continue
            for j in (i+1):length(factors)
                factors[j] isa SUNF || continue
                result = _contract_ff(factors[i]::SUNF, factors[j]::SUNF, N)
                if result !== nothing
                    new_coeff, new_factor = result
                    deleteat!(factors, j)
                    deleteat!(factors, i)
                    coeff_acc = mul_coeff(coeff_acc, new_coeff)
                    new_factor !== nothing && push!(factors, new_factor)
                    changed = true
                    break
                end
            end
            changed && break
        end
    end

    alg_from_factors(factors, normalise_coeff(coeff_acc))
end

# ---- Delta substitution into another factor ----
# δ^{ab} * F(a,...) → F(b,...) (replace a with b in F)
function _delta_substitute(d::SUNDelta, f::SUNF)
    _replace_adj(f, d.a, d.b) !== nothing && return _replace_adj(f, d.a, d.b)
    _replace_adj(f, d.b, d.a)
end
function _delta_substitute(d::SUNDelta, f::SUND)
    _replace_adj_d(f, d.a, d.b) !== nothing && return _replace_adj_d(f, d.a, d.b)
    _replace_adj_d(f, d.b, d.a)
end
function _delta_substitute(d::SUNDelta, g::SUNDelta)
    # δ^{ab} applied to δ^{ac} → δ^{bc}
    # Must share at least one index
    if g.a == d.a
        return SUNDelta(d.b, g.b)
    elseif g.a == d.b
        return SUNDelta(d.a, g.b)
    elseif g.b == d.a
        return SUNDelta(g.a, d.b)
    elseif g.b == d.b
        return SUNDelta(g.a, d.a)
    end
    nothing
end
_delta_substitute(::SUNDelta, ::AlgFactor) = nothing

function _replace_adj(f::SUNF, old::AdjointIndex, new_idx::AdjointIndex)
    f.a == old && return SUNF(new_idx, f.b, f.c)
    f.b == old && return SUNF(f.a, new_idx, f.c)
    f.c == old && return SUNF(f.a, f.b, new_idx)
    nothing
end

function _replace_adj_d(d::SUND, old::AdjointIndex, new_idx::AdjointIndex)
    d.a == old && return SUND(new_idx, d.b, d.c)
    d.b == old && return SUND(d.a, new_idx, d.c)
    d.c == old && return SUND(d.a, d.b, new_idx)
    nothing
end

# ---- f^{acd} f^{bcd} = N δ^{ab} (2 shared indices) ----
function _contract_ff(f1::SUNF, f2::SUNF, N::Int)
    shared, free1, free2 = _shared_adj([f1.a, f1.b, f1.c], [f2.a, f2.b, f2.c])
    length(shared) == 2 || return nothing
    sign = f1.sign * f2.sign
    (Rational{Int}(sign * N), SUNDelta(free1[1], free2[1]))
end

function _shared_adj(as::Vector{AdjointIndex}, bs::Vector{AdjointIndex})
    shared = AdjointIndex[]
    free_a = AdjointIndex[]
    matched_b = falses(length(bs))
    for a in as
        found = false
        for (j, b) in enumerate(bs)
            if !matched_b[j] && a == b
                push!(shared, a)
                matched_b[j] = true
                found = true
                break
            end
        end
        found || push!(free_a, a)
    end
    free_b = [bs[j] for j in eachindex(bs) if !matched_b[j]]
    (shared, free_a, free_b)
end
