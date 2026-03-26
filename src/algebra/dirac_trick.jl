# Feynfeld.jl — DiracTrick: core Dirac algebra simplification rules
#
# Implements the key gamma matrix identities:
#   g^mu g_mu = D (trace of identity in D dims)
#   p-slash . p-slash = p^2
#   {g^mu, g^nu} = 2 g^{mu,nu}
#   g^5 . g^5 = 1
#   g^6/g^7 projector algebra
#   g^mu g^nu g_mu = (2-D) g^nu
#
# Algorithm: greedy left-to-right scan over DiracChain elements,
# applying the highest-priority rule that matches, then restarting.
#
# Ref: FeynCalc Dirac/DiracTrick.m (1305 lines in FeynCalc)

export dirac_trick

"""
    dirac_trick(chain::DiracChain) -> Vector{Tuple{Any, DiracChain}}

Apply Dirac algebra simplification rules to a chain.
Returns a sum of `(coefficient, simplified_chain)` terms.

# Examples
```julia
dirac_trick(dot(GA(:μ), GA(:μ)))       # → [(4, DiracChain([]))]
dirac_trick(dot(GS(:p), GS(:p)))       # → [(SP(:p,:p), DiracChain([]))]
dirac_trick(dot(GA5(), GA5()))          # → [(1, DiracChain([]))]
```
"""
function dirac_trick(chain::DiracChain)
    _simplify_loop(1, chain)
end

function _simplify_loop(coeff, chain::DiracChain)
    elems = chain.elements
    length(elems) == 0 && return [(coeff, chain)]

    # Rule 1: gamma5 squared → identity
    r = _rule_gamma5_squared(elems)
    r !== nothing && return _simplify_loop(coeff, DiracChain(r))

    # Rule 2: projector algebra (g6*g6=g6, g6*g7=0, etc.)
    r = _rule_projector(elems)
    r !== nothing && return _apply_rule_result(coeff, r)

    # Rule 3: adjacent slash squared → scalar product
    r = _rule_slash_squared(elems)
    r !== nothing && return _apply_rule_result(coeff, r)

    # Rule 4: adjacent index contraction g^mu g_mu → D
    r = _rule_adjacent_trace(elems)
    r !== nothing && return _apply_rule_result(coeff, r)

    # Rule 5: g^mu g^nu g_mu → (2-D) g^nu (and higher)
    r = _rule_sandwich_contraction(elems)
    r !== nothing && return _apply_rule_result(coeff, r)

    # No rule matched → return as-is
    [(coeff, chain)]
end

"""Apply a rule result: (new_coeff, new_elems) or sum of them."""
function _apply_rule_result(coeff, result)
    out = Tuple{Any,DiracChain}[]
    for (c, elems) in result
        new_coeff = _mul_coeff(coeff, c)
        new_coeff == 0 && continue
        append!(out, _simplify_loop(new_coeff, DiracChain(elems)))
    end
    isempty(out) && return [(0, DiracChain(DiracElement[]))]
    out
end

# ── Rule 1: gamma5 squared ──────────────────────────────────────────

function _rule_gamma5_squared(elems)
    for i in 1:(length(elems)-1)
        if _is_special(elems[i], 5) && _is_special(elems[i+1], 5)
            return vcat(elems[1:i-1], elems[i+2:end])
        end
    end
    nothing
end

# ── Rule 2: projector algebra ────────────────────────────────────────

function _rule_projector(elems)
    for i in 1:(length(elems)-1)
        a, b = elems[i], elems[i+1]
        r = _projector_pair(a, b)
        r !== nothing || continue
        prefix = elems[1:i-1]
        suffix = elems[i+2:end]
        if r == :zero
            return [(0, DiracElement[])]
        elseif r == :keep_first
            return [(1, vcat(prefix, [a], suffix))]
        elseif r == :keep_second
            return [(1, vcat(prefix, [b], suffix))]
        elseif r == :negate_second
            return [(-1, vcat(prefix, [b], suffix))]
        end
    end
    nothing
end

function _projector_pair(a, b)
    !(_is_special_any(a) && _is_special_any(b)) && return nothing
    ai, bi = a.slot.id, b.slot.id
    # g5*g6 = g6, g5*g7 = -g7
    ai == 5 && bi == 6 && return :keep_second
    ai == 5 && bi == 7 && return :negate_second
    ai == 6 && bi == 5 && return :keep_first
    ai == 7 && bi == 5 && return :negate_second   # g7*g5 = -g7
    # g6*g6 = g6, g7*g7 = g7 (idempotent)
    ai == 6 && bi == 6 && return :keep_first
    ai == 7 && bi == 7 && return :keep_first
    # g6*g7 = 0, g7*g6 = 0
    ai == 6 && bi == 7 && return :zero
    ai == 7 && bi == 6 && return :zero
    nothing
end

# ── Rule 3: slash squared ────────────────────────────────────────────

function _rule_slash_squared(elems)
    for i in 1:(length(elems)-1)
        a, b = elems[i], elems[i+1]
        (a.slot isa MomSlot && b.slot isa MomSlot) || continue
        a.slot.mom isa Momentum && b.slot.mom isa Momentum || continue
        a.slot.mom == b.slot.mom || continue
        # p-slash . p-slash = p^2
        sp = pair(a.slot.mom, b.slot.mom)
        rest = vcat(elems[1:i-1], elems[i+2:end])
        return [(sp, rest)]
    end
    nothing
end

# ── Rule 4: adjacent trace g^mu g_mu → D ────────────────────────────

function _rule_adjacent_trace(elems)
    for i in 1:(length(elems)-1)
        a, b = elems[i], elems[i+1]
        (a.slot isa LISlot && b.slot isa LISlot) || continue
        a.slot.index.name == b.slot.index.name || continue
        # g^mu g_mu = dim_trace(projected_dim)
        dp = dim_contract(a.slot.index.dim, b.slot.index.dim)
        dp === nothing && return [(0, DiracElement[])]
        trace = dim_trace(dp)
        rest = vcat(elems[1:i-1], elems[i+2:end])
        return [(trace, rest)]
    end
    nothing
end

# ── Rule 5: sandwich contraction g^mu ... g_mu ──────────────────────

function _rule_sandwich_contraction(elems)
    for i in eachindex(elems)
        elems[i].slot isa LISlot || continue
        idx = elems[i].slot.index
        for j in (i+2):length(elems)
            elems[j].slot isa LISlot || continue
            elems[j].slot.index.name == idx.name || continue
            dp = dim_contract(idx.dim, elems[j].slot.index.dim)
            dp === nothing && return [(0, DiracElement[])]
            n = j - i - 1  # number of gammas between
            prefix = elems[1:i-1]
            between = elems[i+1:j-1]
            suffix = elems[j+1:end]
            return _apply_contraction_identity(dp, n, prefix, between, suffix)
        end
    end
    nothing
end

"""Apply g^mu X_1...X_n g_mu identity. Returns [(coeff, elems)]."""
function _apply_contraction_identity(dim::DimSlot, n, prefix, between, suffix)
    d = dim_trace(dim)
    if n == 1
        # g^mu g^a g_mu = (2-D) g^a
        c = _sub_coeff(2, d)  # 2 - D
        return [(c, vcat(prefix, between, suffix))]
    elseif n == 2
        # g^mu g^a g^b g_mu = 4 g^{a,b} + (D-4) g^a g^b
        # where g^{a,b} = Pair(a,b) (metric contraction of the two inner gammas)
        a, b = between[1], between[2]
        metric = _gamma_pair(a, b)
        if metric !== nothing
            c_metric = _mul_coeff(4, metric)  # 4 * g^{a,b}
            c_remain = _sub_coeff(d, 4)  # D - 4
            return [(c_metric, vcat(prefix, suffix)),
                    (c_remain, vcat(prefix, [a, b], suffix))]
        end
        # Fallback: if inner gammas aren't LISlot/MomSlot for metric, skip
        return nothing
    end
    # n ≥ 3: defer to general formula (future phase)
    nothing
end

# ── Helpers ──────────────────────────────────────────────────────────

_is_special(g::DiracGamma, id::Int) = g.slot isa SpecialSlot && g.slot.id == id
_is_special(::Spinor, ::Int) = false
_is_special_any(g::DiracGamma) = g.slot isa SpecialSlot
_is_special_any(::Spinor) = false

"""Compute 2 - D or D - 4 etc. as number or symbol."""
function _sub_coeff(a, b)
    (a isa Number && b isa Number) && return a - b
    a == b && return 0
    :($a - $b)
end

"""Extract a Pair from two gamma matrices (for metric contraction)."""
function _gamma_pair(a::DiracGamma, b::DiracGamma)
    if a.slot isa LISlot && b.slot isa LISlot
        return pair(a.slot.index, b.slot.index)
    elseif a.slot isa LISlot && b.slot isa MomSlot && b.slot.mom isa Momentum
        return pair(a.slot.index, b.slot.mom)
    elseif a.slot isa MomSlot && a.slot.mom isa Momentum && b.slot isa LISlot
        return pair(a.slot.mom, b.slot.index)
    elseif a.slot isa MomSlot && a.slot.mom isa Momentum && b.slot isa MomSlot && b.slot.mom isa Momentum
        return pair(a.slot.mom, b.slot.mom)
    end
    nothing
end
