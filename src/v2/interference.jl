# Cross-line traces and interference for multi-channel amplitudes.
#
# Generalizes spin_sum.jl's _single_line_trace to handle:
# - Cross terms between two chains on the same fermion line (Compton)
# - Reconnected fermion lines from different channels (Bhabha s×t)
#
# Uses _completeness, _conjugate_gammas, _build_and_trace from spin_sum.jl.

# ---- Cross-line trace ----

# Compute Tr[(p̸_R+m_R) Γ̃_conj (p̸_L+m_L) Γ_fwd] where Γ_fwd comes from
# chain_fwd and Γ̃_conj comes from chain_conj (reversed + relabeled).
# Both chains must share the same spinor momenta at boundaries.
# Reduces to _single_line_trace when chain_fwd == chain_conj.
function _cross_line_trace(chain_fwd::DiracChain, chain_conj::DiracChain)
    elems_f = chain_fwd.elements
    sp_L = first(elems_f)::Spinor
    sp_R = last(elems_f)::Spinor
    gammas_f = DiracGamma[e for e in elems_f[2:end-1] if e isa DiracGamma]

    elems_c = chain_conj.elements
    gammas_c = DiracGamma[e for e in elems_c[2:end-1] if e isa DiracGamma]
    gammas_c_rev = _conjugate_gammas(gammas_c)

    _, mass_R, mom_R = _completeness(sp_R)
    _, mass_L, mom_L = _completeness(sp_L)

    _build_and_trace(mom_R, mass_R, gammas_c_rev, mom_L, mass_L, gammas_f)
end

# ---- Boson-exchange interference: reconnected fermion lines ----

"""
    spin_sum_interference(amp_i, amp_j)

Spin-summed cross term Σ_spins M_i* M_j for two boson-exchange channels
with reconnected fermion lines (e.g. Bhabha s×t).

Returns a single long trace (AlgSum). For the interference term in
|M|² = ... - 2·Re(M_i* M_j)/(d_i·d_j), use this as T_int.
"""
function spin_sum_interference(
    amp_i::NTuple{2,DiracChain},
    amp_j::NTuple{2,DiracChain},
)
    iL, iR = _line_info(amp_i[1]), _line_info(amp_i[2])
    jL, jR = _line_info(amp_j[1]), _line_info(amp_j[2])

    # For M_i*, the spinor roles reverse (bar <-> plain).
    # Original iL: v̄(p2) Γ u(p1) → Conjugate: ū(p1) Γ̃ v(p2)
    i_conj = [_SpinSumLineInfo(l.plain_mom, l.bar_mom, l.gammas) for l in [iL, iR]]
    j_lines = [jL, jR]

    # Follow the loop: j-sink → i*-gammas → j-source → j-gammas → ...
    trace_gs = DiracGamma[]
    current_mom = jL.plain_mom
    visited = Set{Symbol}()

    while current_mom.name ∉ visited
        push!(visited, current_mom.name)
        push!(trace_gs, GS(current_mom))

        i_line = _find_line_by_bar_mom(i_conj, current_mom)
        append!(trace_gs, _conjugate_gammas(i_line.gammas))

        next_mom = i_line.plain_mom
        push!(trace_gs, GS(next_mom))

        j_line = _find_line_by_bar_mom(j_lines, next_mom)
        append!(trace_gs, j_line.gammas)

        current_mom = j_line.plain_mom
    end

    dirac_trace(trace_gs)
end

# ---- Helpers ----

struct _SpinSumLineInfo
    bar_mom::Momentum
    plain_mom::Momentum
    gammas::Vector{DiracGamma}
end

function _line_info(chain::DiracChain)
    elems = chain.elements
    sp_bar = first(elems)::Spinor
    sp_plain = last(elems)::Spinor
    gs = DiracGamma[e for e in elems[2:end-1] if e isa DiracGamma]
    _SpinSumLineInfo(sp_bar.momentum, sp_plain.momentum, gs)
end

function _find_line_by_bar_mom(lines::Vector{_SpinSumLineInfo}, mom::Momentum)
    idx = findfirst(l -> l.bar_mom.name == mom.name, lines)
    idx === nothing && error("No line with bar momentum $(mom.name)")
    lines[idx]
end
