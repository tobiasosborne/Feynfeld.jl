# Feynfeld.jl — Dirac equation simplification
#
# Applies the Dirac equation at chain boundaries:
#   p-slash . u(p) = m * u(p)
#   ubar(p) . p-slash = m * ubar(p)
#   p-slash . v(p) = -m * v(p)
#   vbar(p) . p-slash = -m * vbar(p)
#
# Ref: FeynCalc Dirac/DiracEquation.m

export dirac_equation

"""
    dirac_equation(chain::DiracChain) -> Vector{Tuple{Any, DiracChain}}

Apply the Dirac equation at chain boundaries:
- `p-slash u(p) = m u(p)` at the right end
- `ūbar(p) p-slash = m ūbar(p)` at the left end

Returns `(coefficient, simplified_chain)` pairs.
"""
function dirac_equation(chain::DiracChain)
    elems = chain.elements
    n = length(elems)
    n < 2 && return [(1, chain)]

    # Check right boundary: ... p-slash u(p) or ... p-slash v(p)
    if elems[end] isa Spinor && n >= 2 && elems[end-1] isa DiracGamma
        g = elems[end-1]
        sp = elems[end]
        r = _dirac_eq_right(g, sp)
        if r !== nothing
            mass_coeff, new_spinor = r
            new_elems = vcat(elems[1:end-2], [new_spinor])
            result = [(mass_coeff, DiracChain(new_elems))]
            # Recurse in case there's another p-slash before
            return _apply_deq_recursive(result)
        end
    end

    # Check left boundary: ubar(p) p-slash ... or vbar(p) p-slash ...
    if elems[1] isa Spinor && n >= 2 && elems[2] isa DiracGamma
        sp = elems[1]
        g = elems[2]
        r = _dirac_eq_left(sp, g)
        if r !== nothing
            mass_coeff, new_spinor = r
            new_elems = vcat([new_spinor], elems[3:end])
            result = [(mass_coeff, DiracChain(new_elems))]
            return _apply_deq_recursive(result)
        end
    end

    [(1, chain)]
end

function _apply_deq_recursive(terms)
    out = Tuple{Any,DiracChain}[]
    for (c, ch) in terms
        for (c2, ch2) in dirac_equation(ch)
            push!(out, (_mul_coeff(c, c2), ch2))
        end
    end
    out
end

"""Apply Dirac equation at right boundary: g . spinor."""
function _dirac_eq_right(g::DiracGamma, sp::Spinor)
    g.slot isa MomSlot || return nothing
    g.slot.mom isa Momentum || return nothing
    g.slot.mom == sp.momentum || return nothing
    # p-slash u(p) = m u(p), p-slash v(p) = -m v(p)
    sign = sp.kind in (:u, :ubar) ? 1 : -1
    (_mul_coeff(sign, sp.mass), sp)
end

"""Apply Dirac equation at left boundary: spinor . g."""
function _dirac_eq_left(sp::Spinor, g::DiracGamma)
    g.slot isa MomSlot || return nothing
    g.slot.mom isa Momentum || return nothing
    g.slot.mom == sp.momentum || return nothing
    # ubar(p) p-slash = m ubar(p), vbar(p) p-slash = -m vbar(p)
    sign = sp.kind in (:u, :ubar) ? 1 : -1
    (_mul_coeff(sign, sp.mass), sp)
end
