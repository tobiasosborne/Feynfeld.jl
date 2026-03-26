# Feynfeld.jl — Non-commutative Dirac chain (DOT)
#
# DiracChain is the non-commutative product of gamma matrices and spinors.
# It is the Julia analog of FeynCalc's DOT[a, b, c, ...].
#
# Key operations:
#   dot(a, b, ...) — construct a chain, flattening nested chains
#   dot_simplify(chain) — expand MomentumSum inside gammas, distribute over sums
#
# Ref: FeynCalc NonCommAlgebra/DotSimplify.m

export dot, dot_simplify, DiracElement

const DiracElement = Union{DiracGamma,Spinor}

"""
    dot(a, b, ...) -> DiracChain

Construct a non-commutative product of DiracGamma/Spinor objects.
Flattens nested DiracChains. Filters out identity (nothing) elements.

# Examples
```julia
dot(GA(:μ), GA(:ν))                  # γ^μ γ^ν
dot(GA(:μ), GS(:p), GA(:ν))         # γ^μ p̸ γ^ν
dot(Spinor(:ubar, Momentum(:p), :m), GA(:μ), Spinor(:u, Momentum(:k), :m))
```
"""
function dot(args::Union{DiracElement,DiracChain}...)
    elements = DiracElement[]
    for a in args
        if a isa DiracChain
            append!(elements, a.elements)
        else
            push!(elements, a)
        end
    end
    DiracChain(elements)
end

# Allow building chains incrementally
Base.:*(a::DiracElement, b::DiracElement) = dot(a, b)
Base.:*(a::DiracChain, b::DiracElement) = dot(a, b)
Base.:*(a::DiracElement, b::DiracChain) = dot(a, b)
Base.:*(a::DiracChain, b::DiracChain) = dot(a, b)

# ── DotSimplify ──────────────────────────────────────────────────────

"""
    dot_simplify(chain::DiracChain) -> Vector{Tuple{Any, DiracChain}}

Expand MomentumSum slots inside DiracGamma elements of a chain,
distributing into a sum of chains.

Returns a list of `(coefficient, chain)` pairs.
A chain with no MomentumSum returns `[(1, chain)]` unchanged.

# Example
```julia
# γ^μ (p+q)-slash → γ^μ p-slash + γ^μ q-slash
pq = Momentum(:p) + Momentum(:q)
chain = dot(GA(:μ), dirac_gamma(pq))
dot_simplify(chain)
# → [(1//1, dot(GA(:μ), GS(:p))), (1//1, dot(GA(:μ), GS(:q)))]
```
"""
function dot_simplify(chain::DiracChain)
    # Find first element with MomentumSum
    idx = _find_momentum_sum(chain)
    idx === nothing && return [(1 // 1, chain)]

    # Expand the MomentumSum element into individual terms
    g = chain.elements[idx]
    ms = g.slot.mom  # ::MomentumSum
    prefix = chain.elements[1:idx-1]
    suffix = chain.elements[idx+1:end]

    result = Tuple{Any,DiracChain}[]
    for (coeff, name) in ms.terms
        expanded_g = DiracGamma(MomSlot(Momentum(name, ms.dim)))
        new_chain = DiracChain(vcat(prefix, [expanded_g], suffix))
        # Recursively expand remaining MomentumSums
        for (c2, ch2) in dot_simplify(new_chain)
            push!(result, (coeff * c2, ch2))
        end
    end

    isempty(result) && return [(0 // 1, DiracChain(DiracElement[]))]
    result
end

"""Find the first DiracGamma with a MomentumSum slot, or nothing."""
function _find_momentum_sum(chain::DiracChain)
    for (i, e) in enumerate(chain.elements)
        if e isa DiracGamma && e.slot isa MomSlot && e.slot.mom isa MomentumSum
            return i
        end
    end
    nothing
end
