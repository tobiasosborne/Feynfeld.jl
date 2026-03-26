# Feynfeld.jl — FermionSpinSum: completeness relations for spin-averaged |M|²
#
# Applies: Σ_s u(p,s) ū(p,s) = p̸ + m
#          Σ_s v(p,s) v̄(p,s) = p̸ - m
#
# Takes two matched spinor chains (from M and M†), inserts completeness
# relations at the junction points, and returns the Dirac trace as AlgSum.
#
# Ref: FeynCalc Feynman/FermionSpinSum.m, Peskin & Schroeder §5.1

export fermion_spin_sum

"""
    fermion_spin_sum(chain_M, chain_Mconj) -> AlgSum

Apply completeness relations to a pair of spinor chains forming a closed
fermion loop (from |M|² spin averaging).

  chain_M:     [ū(k), Γ₁, ..., Γₙ, v(k')]     (from amplitude M)
  chain_Mconj: [v̄(k'), Γ'₁, ..., Γ'ₘ, u(k)]   (from M†)

Returns: Σ over completeness terms of Tr[(k̸+m) Γ₁...Γₙ (k̸'-m') Γ'₁...Γ'ₘ]
as an AlgSum (traces evaluated).

For massless fermions, this is a single trace. For massive, 4 terms.
"""
function fermion_spin_sum(chain_M::DiracChain, chain_Mconj::DiracChain)
    elems_M = chain_M.elements
    elems_Mc = chain_Mconj.elements

    length(elems_M) < 2 && error("chain_M must have at least 2 elements (two spinors)")
    length(elems_Mc) < 2 && error("chain_Mconj must have at least 2 elements (two spinors)")

    sp_L = elems_M[1]      # left spinor of M (ū or v̄)
    sp_R = elems_M[end]    # right spinor of M (u or v)
    sp_Lc = elems_Mc[1]    # left spinor of M† (conjugate of sp_R)
    sp_Rc = elems_Mc[end]  # right spinor of M† (conjugate of sp_L)

    for s in (sp_L, sp_R, sp_Lc, sp_Rc)
        s isa Spinor || error("Chain endpoints must be Spinors, got $(typeof(s))")
    end

    _validate_spinor_matching(sp_R, sp_Lc)
    _validate_spinor_matching(sp_Rc, sp_L)

    inner_M = elems_M[2:end-1]
    inner_Mc = elems_Mc[2:end-1]

    # Completeness at left junction: sp_Rc...sp_L → Σ u(k)ū(k) or v(k)v̄(k)
    comp_L = _completeness(sp_L)
    # Completeness at right junction: sp_R...sp_Lc → same
    comp_R = _completeness(sp_R)

    # Distribute: Tr[comp_L · inner_M · comp_R · inner_Mc]
    result = alg_zero()
    for (c_L, gammas_L) in comp_L
        for (c_R, gammas_R) in comp_R
            chain_elems = Union{DiracGamma,Spinor}[
                gammas_L; inner_M; gammas_R; inner_Mc
            ]
            chain = DiracChain(chain_elems)
            tr = dirac_trace_alg(chain)
            coeff = _mul_coeff(c_L, c_R)
            result = result + _scale(tr, coeff)
        end
    end
    result
end

# ── Completeness relation ──────────────────────────────────────────

"""
Return completeness insertion as [(coeff, DiracGamma[])...].
  u-type: p̸ + m → [(1, [GS(p)]), (m, [])]
  v-type: p̸ - m → [(1, [GS(p)]), (-m, [])]
If mass is 0, returns just [(1, [GS(p)])].
"""
function _completeness(sp::Spinor)
    p_slash = DiracGamma(MomSlot(sp.momentum))
    sign = sp.kind in (:u, :ubar) ? 1 : -1

    terms = Tuple{Any,Vector{DiracGamma}}[(1, [p_slash])]
    if sp.mass != 0
        push!(terms, (_mul_coeff(sign, sp.mass), DiracGamma[]))
    end
    terms
end

# ── Validation ─────────────────────────────────────────────────────

"""Validate that two spinors at a junction are conjugate partners."""
function _validate_spinor_matching(sp_out::Spinor, sp_in::Spinor)
    sp_out.momentum == sp_in.momentum || error(
        "Spinor momentum mismatch: $(sp_out.momentum) vs $(sp_in.momentum)")
    _is_conjugate_pair(sp_out.kind, sp_in.kind) || error(
        "Spinor kinds don't form a conjugate pair: :$(sp_out.kind) and :$(sp_in.kind)")
end

"""Check if two spinor kinds form a valid completeness pair."""
function _is_conjugate_pair(a::Symbol, b::Symbol)
    (a == :u && b == :ubar) || (a == :v && b == :vbar) ||
    (a == :ubar && b == :u) || (a == :vbar && b == :v)
end
