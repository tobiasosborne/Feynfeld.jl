# Feynfeld.jl — Lorentz index contraction
#
# Ports FeynCalc's Contract (Contract.m, PairContract.m).
# Algorithm: worklist reducer that finds repeated LorentzIndex names
# across products of Pairs and performs Einstein summation.
#
# Key contraction rules:
#   g^{μν} g_{νρ} = g^{μ}_{ρ}           (metric-metric)
#   g^{μν} p_ν    = p^μ                  (metric-vector)
#   p^μ q_μ       = p·q                  (vector-vector → scalar product)
#   g^{μ}_{μ}     = dim_trace(dim)       (metric trace)

export contract

"""
    contract(pairs...; ctx=SPContext()) -> result

Contract all repeated Lorentz indices in a product of Pairs.

Returns the simplified result: a scalar (number/symbol), a single Pair,
or a tuple `(coefficient, remaining_pairs)` for partially contracted products.

# Examples
```julia
contract(MT(:μ, :μ))                     # → 4
contract(FV(:p, :μ), FV(:q, :μ))        # → SP(:p, :q)
contract(MT(:μ, :ν), MT(:ν, :ρ))        # → MT(:μ, :ρ)
contract(MTD(:μ, :μ))                    # → :D
```
"""
function contract(p::Feynfeld.Pair; ctx::SPContext=SPContext())
    coeff, remaining = _contract_product([p], ctx)
    _simplify_result(coeff, remaining)
end

function contract(p1::Feynfeld.Pair, p2::Feynfeld.Pair, ps::Feynfeld.Pair...;
                  ctx::SPContext=SPContext())
    all_pairs = Feynfeld.Pair[p1, p2, ps...]
    coeff, remaining = _contract_product(all_pairs, ctx)
    _simplify_result(coeff, remaining)
end

# ── Core algorithm ───────────────────────────────────────────────────

"""Worklist contraction: repeatedly find and contract dummy index pairs."""
function _contract_product(factors::Vector{Feynfeld.Pair}, ctx::SPContext)
    factors = copy(factors)
    coeff::Any = 1

    changed = true
    while changed
        changed = false
        inv = _index_inventory(factors)

        for (name, positions) in inv
            length(positions) < 2 && continue
            length(positions) > 2 && error(
                "Index :$name appears $(length(positions)) times — " *
                "violates Einstein summation convention"
            )

            i, j = positions[1], positions[2]

            if i == j
                # Self-contraction: g^{μ}_{μ} = trace
                coeff = _mul_coeff(coeff, dim_trace(factors[i].a.dim))
                deleteat!(factors, i)
            else
                # Cross-pair contraction
                surv_i = _surviving_arg(factors[i], name)
                surv_j = _surviving_arg(factors[j], name)

                hi, lo = max(i, j), min(i, j)
                deleteat!(factors, hi)
                deleteat!(factors, lo)

                result = pair(surv_i, surv_j)
                if result == 0
                    return (0, Feynfeld.Pair[])
                end

                result = _apply_sp(result, ctx)
                if result isa Feynfeld.Pair
                    push!(factors, result)
                else
                    coeff = _mul_coeff(coeff, result)
                end
            end

            changed = true
            break  # restart scan — indices shifted
        end
    end

    (coeff, factors)
end

# ── Helpers ──────────────────────────────────────────────────────────

"""Build index inventory: map LorentzIndex name → positions in factors."""
function _index_inventory(factors::Vector{Feynfeld.Pair})
    inv = Dict{Symbol,Vector{Int}}()
    for (i, p) in enumerate(factors)
        for arg in (p.a, p.b)
            if arg isa LorentzIndex
                push!(get!(inv, arg.name, Int[]), i)
            end
        end
    end
    inv
end

"""Get the slot that is NOT the contracted index."""
function _surviving_arg(p::Feynfeld.Pair, idx_name::Symbol)
    a_match = p.a isa LorentzIndex && p.a.name == idx_name
    b_match = p.b isa LorentzIndex && p.b.name == idx_name
    a_match && b_match && return nothing  # self-contraction
    a_match && return p.b
    b_match && return p.a
    error("BUG: index :$idx_name not found in Pair")
end

"""Apply SP registry: if Pair is Momentum-Momentum, check for assigned value."""
function _apply_sp(p::Feynfeld.Pair, ctx::SPContext)
    if p.a isa Momentum && p.b isa Momentum
        val = get_sp(ctx, p.a.name, p.b.name)
        val !== nothing && return val
    end
    p
end

"""Multiply coefficients, handling symbolic values (e.g. :D)."""
function _mul_coeff(a, b)
    a == 0 && return 0
    b == 0 && return 0
    a == 1 && return b
    b == 1 && return a
    (a isa Number && b isa Number) && return a * b
    :($a * $b)
end

"""Simplify the (coeff, pairs) tuple into a natural return value."""
function _simplify_result(coeff, remaining::Vector{Feynfeld.Pair})
    if isempty(remaining)
        return coeff
    elseif coeff == 1 && length(remaining) == 1
        return remaining[1]
    elseif coeff == 0
        return 0
    elseif coeff == 1
        return Tuple(remaining)
    else
        return (coeff, Tuple(remaining))
    end
end

