# Feynfeld.jl — SU(N) colour trace evaluation
#
# Tr(1) = N
# Tr(T^a) = 0
# Tr(T^a T^b) = (1/2) δ^{ab}
# Tr(T^a T^b T^c) = (1/4)(d^{abc} + i f^{abc})
# Recursive: Tr(T^a T^b rest) = δ^{ab}/(2N) Tr(rest) + (1/2)d^{abf} Tr(T^f rest) + (i/2)f^{abf} Tr(T^f rest)
#
# Ref: FeynCalc SUN/SUNTrace.m

export SUNTrace, sun_trace

"""Unevaluated colour trace wrapper."""
struct SUNTrace <: FeynExpr
    chain::ColourChain
end

# Counter for fresh dummy indices
const _SUN_DUMMY_COUNTER = Ref(0)
function _fresh_sun_index()
    _SUN_DUMMY_COUNTER[] += 1
    SUNIndex(Symbol("_sun_", _SUN_DUMMY_COUNTER[]))
end

"""
    sun_trace(generators::SUNT...) -> result
    sun_trace(chain::ColourChain) -> result

Evaluate the colour trace of a product of SU(N) generators.

# Examples
```julia
sun_trace()                          # Tr(1) = :N
sun_trace(SUNT(SUNIndex(:a)))        # Tr(T^a) = 0
sun_trace(SUNT(SUNIndex(:a)), SUNT(SUNIndex(:b)))  # (1//2) SUNDelta
```
"""
function sun_trace(chain::ColourChain)
    _sun_trace_eval(chain.elements)
end

function sun_trace(gs::SUNT...)
    _sun_trace_eval(collect(gs))
end

function _sun_trace_eval(elems::Vector{SUNT})
    n = length(elems)

    # Tr(1) = N
    n == 0 && return SUNN

    # Tr(T^a) = 0
    n == 1 && return 0

    # Tr(T^a T^b) = (1//2) δ^{ab}
    if n == 2
        return (1 // 2, SUNDelta(elems[1].a, elems[2].a))
    end

    # Tr(T^a T^b T^c) = (1//4)(d^{abc} + i*f^{abc})
    if n == 3
        a, b, c = elems[1].a, elems[2].a, elems[3].a
        return (1 // 4, SUND(a, b, c), SUNF(a, b, c))
    end

    # Recursive: peel first two generators
    # T^a T^b = δ^{ab}/(2N) + (1/2)(d^{abf} + i*f^{abf}) T^f
    a = elems[1].a
    b = elems[2].a
    rest = elems[3:end]
    f = _fresh_sun_index()

    # Term 1: δ^{ab}/(2N) * Tr(rest)
    sub1 = _sun_trace_eval(rest)

    # Term 2: (1/2) d^{abf} * Tr(T^f rest)
    sub2 = _sun_trace_eval(vcat([SUNT(f)], rest))

    # Term 3: (i/2) f^{abf} * Tr(T^f rest)
    # (same sub-trace as term 2, different structure constant)
    sub3 = sub2  # reuse

    return (SUNDelta(a, b), sub1, SUND(a, b, f), SUNF(a, b, f), sub2)
end
