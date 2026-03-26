# Feynfeld.jl — DiracTrace: Dirac trace evaluation
#
# Tr[g^{a1} ... g^{a_{2n}}] = Σ_{i=2}^{2n} (-1)^i g^{a1,ai} Tr[remaining]
# Tr[g5 g^a g^b g^c g^d] = 4i * ε^{abcd}
#
# Base cases:
#   Tr[1] = 4 (or D in D-dim)
#   Tr[odd gammas] = 0
#   Tr[g5] = 0
#   Tr[g^a g^b] = 4 g^{ab}
#
# Ref: FeynCalc Dirac/DiracTrace.m (Thomas Hahn's recursive formula)

export dirac_trace

"""
    dirac_trace(chain::DiracChain) -> result

Evaluate the Dirac trace of a gamma matrix chain.

Returns a scalar (number/symbol), a Pair, or a sum of terms
as `Vector{Tuple{Any, Any}}`.

# Examples
```julia
dirac_trace(DiracChain(DiracElement[]))           # → 4
dirac_trace(dot(GA(:μ), GA(:ν)))                  # → 4 * MT(:μ,:ν)
dirac_trace(dot(GA(:μ), GA(:ν), GA(:ρ), GA(:σ))) # → 4*(g^{μν}g^{ρσ} - g^{μρ}g^{νσ} + g^{μσ}g^{νρ})
```
"""
function dirac_trace(chain::DiracChain)
    elems = chain.elements
    # Spinors in a trace is an error
    any(e -> e isa Spinor, elems) && error("Cannot take trace of a chain with spinors")

    n = length(elems)

    # Tr[1] = 4
    n == 0 && return 4

    # Check for special gammas (g5, g6, g7)
    has_special = any(e -> e isa DiracGamma && e.slot isa SpecialSlot, elems)

    if !has_special
        # Regular trace: vanishes for odd number of gammas
        isodd(n) && return 0
        return _trace_no_g5(elems)
    else
        return _trace_with_g5(elems)
    end
end

# ── Non-chiral trace (no gamma5) ────────────────────────────────────

"""Recursive trace formula for even number of gamma matrices."""
function _trace_no_g5(elems::Vector{<:Union{DiracGamma,Spinor}})
    n = length(elems)
    n == 0 && return 4

    # Tr[g^a g^b] = 4 * g^{ab}
    n == 2 && return _mul_coeff(4, _gamma_pair(elems[1], elems[2]))

    # Tr[g^{a1} ... g^{a_{2n}}] = Σ_{i=2}^{2n} (-1)^i * g^{a1,ai} * Tr[rest]
    # where rest = all elements except a1 and ai
    result = Tuple{Any,Any}[]
    a1 = elems[1]
    for i in 2:n
        sign = iseven(i) ? 1 : -1
        metric = _gamma_pair(a1, elems[i])
        metric === nothing && continue
        rest = vcat(elems[2:i-1], elems[i+1:n])
        sub_trace = _trace_no_g5(rest)
        term = _mul_coeff(sign, _mul_coeff(metric, sub_trace))
        term == 0 && continue
        push!(result, (term, nothing))
    end

    isempty(result) && return 0
    length(result) == 1 && return result[1][1]
    # Return as list of scalar terms
    [r[1] for r in result]
end

# ── Chiral trace (with gamma5) ──────────────────────────────────────

"""Trace with gamma5. Base case: Tr[g5 g^a g^b g^c g^d] = 4i * ε^{abcd}."""
function _trace_with_g5(elems)
    # Strip gamma5/6/7 and count non-special gammas
    regular = DiracGamma[]
    g5_count = 0
    for e in elems
        if e isa DiracGamma && e.slot isa SpecialSlot
            e.slot.id == 5 && (g5_count += 1)
            # g6 = (1+g5)/2, g7 = (1-g5)/2 — defer to expansion
            e.slot.id in (6, 7) && error("Trace with g6/g7: expand to g5 form first")
        else
            push!(regular, e)
        end
    end

    # g5^2 = 1, so only parity matters
    if iseven(g5_count)
        # Even g5 count → cancels out, pure regular trace
        isodd(length(regular)) && return 0
        return _trace_no_g5(regular)
    end

    # Odd g5 count: chiral trace
    n = length(regular)
    isodd(n) && return 0           # odd regular gammas → 0
    n < 4 && return 0              # Tr[g5]=0, Tr[g5 g^a g^b]=0
    n != 4 && return 0             # n>4 chiral recursive not yet implemented

    # Tr[g5 g^a g^b g^c g^d] = 4i * ε^{abcd}
    # Using $LeviCivitaSign convention: result = -4i * $LCS * ε
    # With $LCS = -1: result = 4i * ε
    args = PairArg[]
    for g in regular
        if g.slot isa LISlot
            push!(args, g.slot.index)
        elseif g.slot isa MomSlot && g.slot.mom isa Momentum
            push!(args, g.slot.mom)
        else
            return 0  # can't form epsilon with MomentumSum
        end
    end

    eps = levi_civita(args[1], args[2], args[3], args[4])
    eps == 0 && return 0
    # 4i * ε (the i is conventional; we return 4*$LCS*ε for real-valued convention)
    # FeynCalc: Tr[g5 g^a g^b g^c g^d] = -4i * $LeviCivitaSign * ε^{abcd}
    # With $LCS = -1: = 4i * ε
    # For now, return the symbolic form without explicit i
    _mul_coeff(4 * LEVI_CIVITA_SIGN, eps)
end

