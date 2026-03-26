# Feynfeld.jl — PaVeReduce: reduce PaVe coefficients to scalar integrals
#
# Key B-function reductions (Denner 1993):
#   B1  = (1/2p²)(A0(m1²) - A0(m0²) - (p²+m0²-m1²) B0)
#   B00 = (1/6)(A0(m1²) + 2m0² B0 + (p²+m0²-m1²) B1 + m0²+m1²-p²/3)
#   B11 = (1/3p²)((m1²-m0²+p²) B1 - A0(m1²)/2 + A0(m0²)/2 + 2 B00)
#
# Phase 1d: B-function reductions only. C/D reductions deferred to Phase 2.
#
# Ref: FeynCalc LoopIntegrals/PaVeReduce.m, Denner (1993) Eqs (4.18)-(4.20)

export pave_reduce

"""
    pave_reduce(p::PaVe{N}) -> result

Reduce a PaVe coefficient function to scalar integrals where possible.
Currently implements B-function reductions (N=2).

Returns the PaVe unchanged if no reduction is available.
"""
function pave_reduce(p::PaVe{N}) where {N}
    N == 2 && return _reduce_b(p)
    p  # no reduction available for other N yet
end

"""Reduce B-function coefficients to A0 and B0."""
function _reduce_b(p::PaVe{2})
    pp = p.invariants[1]
    m0 = p.masses[1]
    m1 = p.masses[2]
    indices = p.indices

    # B0: already scalar, no reduction
    isempty(indices) && return p

    # B1: expressed in terms of A0, B0
    if indices == [1]
        a0_m0 = A0(m0)
        a0_m1 = A0(m1)
        b0 = B0(pp, m0, m1)
        return (:B1_reduced, a0_m0, a0_m1, b0, pp, m0, m1)
    end

    # B00: expressed in terms of A0, B0, B1
    if indices == [0, 0]
        return (:B00_reduced, A0(m0), A0(m1), B0(pp, m0, m1), pp, m0, m1)
    end

    # B11: expressed in terms of A0, B00, B1
    if indices == [1, 1]
        return (:B11_reduced, A0(m0), A0(m1), B0(pp, m0, m1), pp, m0, m1)
    end

    p  # unrecognized index pattern
end
