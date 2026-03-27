# SPContext: scalar product value storage.
# Uses ScopedValues (Julia 1.11+) for implicit threading + explicit override.
# Values are Coeff (Rational{Int} or DimPoly), not Any.

using Base.ScopedValues

struct SPContext
    values::Dict{Tuple{Symbol,Symbol}, Coeff}
end
SPContext() = SPContext(Dict{Tuple{Symbol,Symbol}, Coeff}())

_sp_key(a::Symbol, b::Symbol) = a <= b ? (a, b) : (b, a)

function set_sp(ctx::SPContext, a::Symbol, b::Symbol, val)
    new_vals = copy(ctx.values)
    new_vals[_sp_key(a, b)] = normalise_coeff(val)
    SPContext(new_vals)
end

function get_sp(ctx::SPContext, a::Symbol, b::Symbol)
    get(ctx.values, _sp_key(a, b), nothing)
end

# ---- ScopedValue for implicit context ----
const CURRENT_SP = ScopedValue(SPContext())

function with_sp(f, ctx::SPContext)
    @with CURRENT_SP => ctx f()
end

# Convenience: build context from Base.Pair assignments
function sp_context(assignments::Vararg{Base.Pair})
    ctx = SPContext()
    for (k, v) in assignments
        a, b = k
        ctx = set_sp(ctx, a, b, v)
    end
    ctx
end

# ---- Evaluate scalar products in an AlgSum ----
function evaluate_sp(s::AlgSum; ctx::SPContext=CURRENT_SP[])
    result = Dict{FactorKey, Coeff}()
    for (fk, c) in s.terms
        new_factors = AlgFactor[]
        new_coeff = c
        for f in fk.factors
            val = _try_sp(f, ctx)
            if val !== nothing
                new_coeff = mul_coeff(new_coeff, val)
            else
                push!(new_factors, f)
            end
        end
        new_coeff = normalise_coeff(new_coeff)
        _coeff_iszero(new_coeff) && continue
        nfk = FactorKey(new_factors)
        existing = get(result, nfk, 0//1)
        new_c = normalise_coeff(add_coeff(existing, new_coeff))
        if _coeff_iszero(new_c)
            delete!(result, nfk)
        else
            result[nfk] = new_c
        end
    end
    AlgSum(result)
end

# Dispatch on Pair type — only ScalarProduct has SP values
_try_sp(p::Pair{Momentum, Momentum}, ctx::SPContext) = get_sp(ctx, p.a.name, p.b.name)
_try_sp(::AlgFactor, ::SPContext) = nothing
