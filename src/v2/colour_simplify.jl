# Colour simplification: delta contraction, structure constant identities.
# All operations return AlgSum or modify AlgSum in place.

# ---- Contract colour deltas in an AlgSum ----
# δ^{ab} contracts with any factor sharing index a or b, replacing it.
function contract_colour(s::AlgSum; N::Int=3)
    result = AlgSum()
    for (fk, c) in s.terms
        add!(result, _contract_colour_term(fk.factors, c, N))
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
            sc = _self_contract(f, N)
            sc === nothing && continue
            coeff_acc = mul_coeff(coeff_acc, sc)
            deleteat!(factors, i)
            changed = true
            break
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

        # Structure constant contractions: f·f, d·d, f·d
        for i in 1:length(factors)
            for j in (i+1):length(factors)
                result = _try_struct_contract(factors[i], factors[j], N)
                result === nothing && continue
                new_coeff, new_factor = result
                iszero(new_coeff) && return AlgSum()  # f·d = 0 kills term
                deleteat!(factors, j); deleteat!(factors, i)
                coeff_acc = mul_coeff(coeff_acc, new_coeff)
                new_factor !== nothing && push!(factors, new_factor)
                changed = true
                break
            end
            changed && break
        end
    end

    alg_from_factors(factors, normalise_coeff(coeff_acc))
end

# ---- Self-contraction dispatch (δ^{aa} → N²-1 or N) ----
_self_contract(d::SUNDelta, N::Int) = d.a == d.b ? Rational{Int}(N^2 - 1) : nothing
_self_contract(d::FundDelta, N::Int) = d.i == d.j ? Rational{Int}(N) : nothing
_self_contract(::AlgFactor, ::Int) = nothing

# ---- Delta substitution into another factor ----
# δ^{ab} * F(a,...) → F(b,...) (replace a with b in F)
function _delta_substitute(d::SUNDelta, f::SUNF)
    result = _replace_adj(f, d.a, d.b)
    result !== nothing && return result
    _replace_adj(f, d.b, d.a)
end
function _delta_substitute(d::SUNDelta, f::SUND)
    result = _replace_adj_d(f, d.a, d.b)
    result !== nothing && return result
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

# ---- Structure constant contraction dispatch ----
_try_struct_contract(a::SUNF, b::SUNF, N) = _contract_ff(a, b, N)
_try_struct_contract(a::SUND, b::SUND, N) = _contract_dd(a, b, N)
_try_struct_contract(a::SUNF, b::SUND, N) = _contract_fd(a, b, N)
_try_struct_contract(a::SUND, b::SUNF, N) = _contract_fd(b, a, N)
_try_struct_contract(::AlgFactor, ::AlgFactor, _) = nothing

# ---- f·f contraction ----
# Ref: refs/papers/Shtabovenko2016_FeynCalc9_1601.01167.pdf, Section 3.4
# Cross-check: refs/FeynCalc/Tests/SUN/SUNSimplify.test
# "f^{acd} f^{bcd} = N δ^{ab}", "f^{abc} f^{abc} = N(N²-1)"
function _contract_ff(f1::SUNF, f2::SUNF, N::Int)
    shared, free1, free2 = _shared_adj([f1.a, f1.b, f1.c], [f2.a, f2.b, f2.c])
    sign = f1.sign * f2.sign
    if length(shared) == 3
        return (Rational{Int}(sign * N * (N^2 - 1)), nothing)
    elseif length(shared) == 2
        return (Rational{Int}(sign * N), SUNDelta(free1[1], free2[1]))
    end
    nothing
end

# ---- d·d contraction ----
# Ref: refs/papers/Shtabovenko2016_FeynCalc9_1601.01167.pdf, Section 3.4
# "d^{acd} d^{bcd} = (N²-4)/N δ^{ab}", "d^{abc} d^{abc} = (N²-1)(N²-4)/N"
function _contract_dd(d1::SUND, d2::SUND, N::Int)
    shared, free1, free2 = _shared_adj([d1.a, d1.b, d1.c], [d2.a, d2.b, d2.c])
    if length(shared) == 3
        return (Rational{Int}((N^2 - 1) * (N^2 - 4)) // N, nothing)
    elseif length(shared) == 2
        return (Rational{Int}(N^2 - 4) // N, SUNDelta(free1[1], free2[1]))
    end
    nothing
end

# ---- f·d contraction ----
# Ref: antisymmetric × symmetric contraction vanishes by index symmetry
# "f^{acd} d^{bcd} = 0", "f^{abc} d^{abc} = 0"
function _contract_fd(f1::SUNF, d1::SUND, N::Int)
    shared, _, _ = _shared_adj([f1.a, f1.b, f1.c], [d1.a, d1.b, d1.c])
    length(shared) >= 2 || return nothing
    (0//1, nothing)  # always zero
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
