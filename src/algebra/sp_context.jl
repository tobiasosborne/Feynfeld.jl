# Feynfeld.jl — Scalar product context
#
# Stores assigned scalar product values (e.g., p·p = m²).
# Explicit context object instead of FeynCalc's global DownValues.
# Thread-safe and test-isolated by design.
#
# Ref: FeynCalc ScalarProduct.m, FCClearScalarProducts.m

export SPContext, set_sp, get_sp

"""
    SPContext()

A copy-on-write collection of scalar product assignments.
Each `set_sp` call returns a new context, leaving the original unchanged.

# Examples
```julia
ctx = SPContext()
ctx = set_sp(ctx, :p, :p, :m²)
ctx = set_sp(ctx, :p, :q, :s)
get_sp(ctx, :p, :p)  # => :m²
get_sp(ctx, :q, :p)  # => :s  (symmetric)
get_sp(ctx, :k, :k)  # => nothing
```
"""
struct SPContext
    values::Dict{Tuple{Symbol,Symbol},Any}
end

SPContext() = SPContext(Dict{Tuple{Symbol,Symbol},Any}())

"""
    set_sp(ctx::SPContext, p::Symbol, q::Symbol, val) -> SPContext

Return a new context with the scalar product p·q = val added.
Keys are stored in canonical (sorted) order.
"""
function set_sp(ctx::SPContext, p::Symbol, q::Symbol, val)
    key = p <= q ? (p, q) : (q, p)
    new_vals = copy(ctx.values)
    new_vals[key] = val
    SPContext(new_vals)
end

"""
    get_sp(ctx::SPContext, p::Symbol, q::Symbol) -> Union{Any, Nothing}

Look up a scalar product value. Returns `nothing` if not set.
"""
function get_sp(ctx::SPContext, p::Symbol, q::Symbol)
    key = p <= q ? (p, q) : (q, p)
    get(ctx.values, key, nothing)
end
