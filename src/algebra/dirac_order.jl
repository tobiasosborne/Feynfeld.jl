# Feynfeld.jl — DiracOrder: normal ordering of gamma matrices
#
# Sorts gamma matrices in a DiracChain into canonical order using
# the anticommutation relation: {γ^μ, γ^ν} = 2g^{μν}
#
# γ^ν γ^μ = 2g^{μν} - γ^μ γ^ν  (for μ < ν, swap with sign + metric term)
#
# Ref: FeynCalc Dirac/DiracOrder.m

export dirac_order

"""
    dirac_order(chain::DiracChain) -> Vector{Tuple{Any, DiracChain}}

Sort gamma matrices into canonical (alphabetical by index name) order
using the anticommutation relation {γ^μ, γ^ν} = 2g^{μν}.

Each swap produces a metric tensor term. Returns a sum of
`(coefficient, chain)` pairs.

# Examples
```julia
dirac_order(dot(GA(:ν), GA(:μ)))
# → [(-1, dot(GA(:μ), GA(:ν))), (2*MT(:μ,:ν), DiracChain([]))]
```
"""
function dirac_order(chain::DiracChain)
    _order_pass(1, chain)
end

function _order_pass(coeff, chain::DiracChain)
    elems = chain.elements
    n = length(elems)
    n <= 1 && return [(coeff, chain)]

    # Find first out-of-order adjacent pair (LISlot gammas only)
    for i in 1:(n-1)
        a, b = elems[i], elems[i+1]
        (a isa DiracGamma && a.slot isa LISlot &&
         b isa DiracGamma && b.slot isa LISlot) || continue
        isless(a.slot.index, b.slot.index) && continue
        a.slot.index == b.slot.index && continue  # same index — leave for DiracTrick

        # Out of order: apply {g^a, g^b} = 2g^{ab}
        # g^a g^b = 2g^{ab} - g^b g^a
        metric = _gamma_pair(a, b)
        swapped = vcat(elems[1:i-1], [b, a], elems[i+2:end])
        prefix = vcat(elems[1:i-1], elems[i+2:end])

        result = Tuple{Any,DiracChain}[]
        # Term 1: -1 * swapped chain (continue ordering)
        append!(result, _order_pass(_mul_coeff(-1, coeff), DiracChain(swapped)))
        # Term 2: 2*metric * remaining chain (skip both gammas)
        if metric !== nothing && metric != 0
            metric_coeff = _mul_coeff(coeff, _mul_coeff(2, metric))
            if metric_coeff != 0
                push!(result, (metric_coeff, DiracChain(prefix)))
            end
        end
        return result
    end

    # Already ordered
    [(coeff, chain)]
end
