# Feynfeld.jl — DiracSimplify: master Dirac algebra simplifier
#
# Orchestrates: DotSimplify → DiracTrick → DiracEquation
# This is the top-level entry point for simplifying Dirac expressions.
#
# Ref: FeynCalc Dirac/DiracSimplify.m

export dirac_simplify

"""
    dirac_simplify(chain::DiracChain; ctx=SPContext()) -> Vector{Tuple{Any, DiracChain}}

Master simplification of a Dirac chain. Applies in order:
1. `dot_simplify` — expand MomentumSum
2. `dirac_trick` — core Dirac algebra rules
3. `dirac_equation` — Dirac equation at chain boundaries

Returns a sum of `(coefficient, chain)` pairs.
"""
function dirac_simplify(chain::DiracChain; ctx::SPContext=SPContext())
    # Step 1: Expand MomentumSum inside gammas
    expanded = dot_simplify(chain)

    # Step 2: Apply DiracTrick to each term
    tricked = Tuple{Any,DiracChain}[]
    for (c1, ch1) in expanded
        for (c2, ch2) in dirac_trick(ch1)
            push!(tricked, (_mul_coeff(c1, c2), ch2))
        end
    end

    # Step 3: Apply Dirac equation at boundaries
    result = Tuple{Any,DiracChain}[]
    for (c, ch) in tricked
        for (c2, ch2) in dirac_equation(ch)
            push!(result, (_mul_coeff(c, c2), ch2))
        end
    end

    # Filter zero terms
    filter!(r -> r[1] != 0, result)
    isempty(result) && return [(0, DiracChain(DiracElement[]))]
    result
end
