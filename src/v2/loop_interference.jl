# Tree × loop interference: spin-summed Born-virtual term.
#
# Computes Σ_spins M_tree* × M_box for 2→2 boson-exchange tree amplitudes
# crossed with 1-loop box amplitudes. Both fermion lines remain independent
# (no reconnection), so the result is a product of two traces.
#
# Each trace combines:
#   - Forward gammas from the BOX amplitude (including loop momentum q)
#   - Conjugate gammas from the TREE amplitude (reversed + relabeled indices)
#
# After trace + contract + expand_sp, the result is an AlgSum containing
# SP(q, p_i) and SP(q, q) factors that feed into the TID.
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 7
# Ref: Peskin & Schroeder, Eq. (5.4) for spin sum via completeness relations

"""
    spin_sum_tree_loop_interference(tree_amp, loop_amp)

Spin-summed interference: Σ_spins M_tree* × M_box.

`tree_amp` and `loop_amp` are each `(chain_e::DiracExpr, chain_mu::DiracExpr)`.
Both must be single-term DiracExpr (QED, no chiral structure).

Returns an AlgSum containing loop-momentum-dependent scalar products
SP(q, p_i) and SP(q, q) as factors in the FactorKey.
"""
function spin_sum_tree_loop_interference(
    tree_amp::NTuple{2, DiracExpr},
    loop_amp::NTuple{2, DiracExpr},
)
    tree_e, tree_mu = tree_amp
    loop_e, loop_mu = loop_amp
    # Each fermion line contributes a trace:
    #   Tr[(completeness_R) × conjugate_tree_gammas × (completeness_L) × loop_gammas]
    tr_e = _interference_line_trace(tree_e, loop_e)
    tr_mu = _interference_line_trace(tree_mu, loop_mu)
    tr_e * tr_mu
end

# Single fermion line trace for tree×loop interference.
#
# For single-term DiracExpr (QED), this is equivalent to
# _cross_line_trace(loop_chain, tree_chain) from interference.jl,
# but uses the forward/conjugate split explicitly.
#
# Trace structure:
#   Tr[(p̸_R + m_R) × Γ̃_tree × (p̸_L + m_L) × Γ_loop]
# where Γ̃_tree is tree gammas reversed and index-relabeled.
function _interference_line_trace(tree_de::DiracExpr, loop_de::DiracExpr)
    # Single-term for QED
    length(tree_de.terms) == 1 || error(
        "Multi-term DiracExpr not yet supported in tree×loop interference")
    length(loop_de.terms) == 1 || error(
        "Multi-term DiracExpr not yet supported in tree×loop interference")

    _, tree_chain = tree_de.terms[1]
    _, loop_chain = loop_de.terms[1]

    # Use _cross_line_trace from interference.jl:
    # forward gammas from loop, conjugate gammas from tree.
    _cross_line_trace(loop_chain, tree_chain)
end
