# Feynfeld.jl — AlgSum operations: contract, expand, trace
#
# Extension methods that integrate AlgSum with existing algebra operations.
# Included after alg_expr.jl and contract.jl so all types and base
# functions (_mul_coeff, _contract_product, expand_scalar_product) exist.

export dirac_trace_alg, evaluate_sp

# ── contract(::AlgSum) ─────────────────────────────────────────────

"""
    contract(s::AlgSum; ctx=SPContext()) -> AlgSum

Contract all repeated Lorentz indices within each term of an AlgSum.
Note: Eps factors are passed through — Pair-Eps index contraction is not
yet implemented (tracked as a known limitation).
"""
function contract(s::AlgSum; ctx::SPContext=SPContext())
    result_terms = AlgTerm[]
    for t in s.terms
        t.coeff == 0 && continue
        contracted = _contract_algterm(t, ctx)
        append!(result_terms, contracted)
    end
    _collect_terms(AlgSum(result_terms))
end

"""Contract a single AlgTerm: extract Pairs, contract, rebuild."""
function _contract_algterm(t::AlgTerm, ctx::SPContext)
    pairs = Feynfeld.Pair[]
    others = AlgFactor[]  # Eps factors pass through (no Pair-Eps contraction yet)
    for f in t.factors
        if f isa Feynfeld.Pair
            push!(pairs, f)
        else
            push!(others, f)
        end
    end

    isempty(pairs) && return [t]

    coeff, remaining = _contract_product(pairs, ctx)
    new_coeff = _mul_coeff(t.coeff, coeff)
    new_coeff == 0 && return AlgTerm[]

    factors = AlgFactor[others; remaining]
    [AlgTerm(new_coeff, factors)]
end

# ── expand_scalar_product(::AlgSum) ────────────────────────────────

"""
    expand_scalar_product(s::AlgSum; ctx=SPContext()) -> AlgSum

Expand all MomentumSum-containing Pairs in each term of an AlgSum.
"""
function expand_scalar_product(s::AlgSum; ctx::SPContext=SPContext())
    result_terms = AlgTerm[]
    for t in s.terms
        t.coeff == 0 && continue
        expanded = _expand_algterm(t, ctx)
        append!(result_terms, expanded)
    end
    _collect_terms(AlgSum(result_terms))
end

"""Expand a single AlgTerm: expand each Pair factor, distribute."""
function _expand_algterm(t::AlgTerm, ctx::SPContext)
    others = AlgFactor[f for f in t.factors if !(f isa Feynfeld.Pair)]
    pairs = Feynfeld.Pair[f for f in t.factors if f isa Feynfeld.Pair]

    isempty(pairs) && return [t]

    current = [AlgTerm(t.coeff, copy(others))]
    for p in pairs
        expanded = expand_scalar_product(p; ctx=ctx)
        next = AlgTerm[]
        for base in current
            for (c, val) in expanded
                rc = c isa Rational && isinteger(c) ? Int(c) : c
                new_coeff = _mul_coeff(base.coeff, rc)
                new_coeff == 0 && continue
                new_factors = copy(base.factors)
                if val isa Feynfeld.Pair
                    push!(new_factors, val)
                elseif val != 1 && val != 0
                    new_coeff = _mul_coeff(new_coeff, val)
                end
                push!(next, AlgTerm(new_coeff, new_factors))
            end
        end
        current = next
    end
    current
end

# ── evaluate_sp: substitute scalar products ────────────────────────

"""
    evaluate_sp(s::AlgSum; ctx::SPContext) -> AlgSum

Substitute all Momentum×Momentum Pair factors using SPContext values.
Factors that have SPContext values become part of the coefficient;
those without are kept as factors.
"""
function evaluate_sp(s::AlgSum; ctx::SPContext=SPContext())
    result_terms = AlgTerm[]
    for t in s.terms
        t.coeff == 0 && continue
        new_coeff = t.coeff
        new_factors = AlgFactor[]
        for f in t.factors
            val = _try_sp_lookup(f, ctx)
            if val !== nothing
                new_coeff = _mul_coeff(new_coeff, val)
            else
                push!(new_factors, f)
            end
        end
        new_coeff == 0 && continue
        push!(result_terms, AlgTerm(new_coeff, new_factors))
    end
    _collect_terms(AlgSum(result_terms))
end

"""Try to look up a Pair's value in SPContext. Returns value or nothing."""
function _try_sp_lookup(f::Feynfeld.Pair, ctx::SPContext)
    f.a isa Momentum && f.b isa Momentum || return nothing
    get_sp(ctx, f.a.name, f.b.name)
end

_try_sp_lookup(::Eps, ::SPContext) = nothing

# ── dirac_trace_alg ────────────────────────────────────────────────

"""
    dirac_trace_alg(chain::DiracChain) -> AlgSum

Evaluate Dirac trace, returning result as an AlgSum.
Wraps `dirac_trace()` and normalizes all return forms.
"""
function dirac_trace_alg(chain::DiracChain)
    raw = dirac_trace(chain)
    _to_algsum(raw)
end

_to_algsum(x::Number) = x == 0 ? alg_zero() : alg_scalar(x)
_to_algsum(x::Symbol) = alg_scalar(x)
_to_algsum(x::Feynfeld.Pair) = alg(x)
_to_algsum(x::Eps) = alg(x)

function _to_algsum(x::Expr)
    coeff, factors = _flatten_product(x)
    isempty(factors) ? alg_scalar(coeff) : AlgSum([AlgTerm(coeff, factors)])
end

function _to_algsum(x::Vector)
    isempty(x) && return alg_zero()
    result = _to_algsum(x[1])
    for i in 2:length(x)
        result = result + _to_algsum(x[i])
    end
    result
end

# ── Product decomposition ─────────────────────────────────────────

"""Recursively decompose a _mul_coeff expression into (scalar_coeff, [AlgFactor...])."""
function _flatten_product(x)
    coeff = Ref{Any}(1)
    factors = AlgFactor[]
    _collect_product!(coeff, factors, x)
    (coeff[], factors)
end

_collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Number) = (c[] = _mul_coeff(c[], x))
_collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Symbol) = (c[] = _mul_coeff(c[], x))
_collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Feynfeld.Pair) = push!(fs, x)
_collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Eps) = push!(fs, x)

function _collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Expr)
    if x.head == :call && length(x.args) >= 3 && x.args[1] == :*
        for arg in @view x.args[2:end]
            _collect_product!(c, fs, arg)
        end
    else
        c[] = _mul_coeff(c[], x)
    end
end

function _collect_product!(c::Ref, fs::Vector{AlgFactor}, x::Vector)
    error("Cannot flatten a Vector inside a product expression. " *
          "Use _to_algsum for sum-level decomposition.")
end
