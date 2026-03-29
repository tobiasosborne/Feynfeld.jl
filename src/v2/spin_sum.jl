# Fermion spin sum: completeness relations via dispatch on SpinorKind.
# sum_s u(p,m) ubar(p,m) = p-slash + m
# sum_s v(p,m) vbar(p,m) = p-slash - m

# Completeness relation: dispatch on spinor kind (not symbol checking!)
_completeness(sp::Spinor{UKind}) = (1, sp.mass, sp.momentum)    # p-slash + m
_completeness(sp::Spinor{VKind}) = (1, -sp.mass, sp.momentum)   # p-slash - m
_completeness(sp::Spinor{UBarKind}) = _completeness_bar(sp)
_completeness(sp::Spinor{VBarKind}) = _completeness_bar(sp)

# Bar spinors delegate to their unbarred partner
_completeness_bar(sp::Spinor{UBarKind}) = (1, sp.mass, sp.momentum)
_completeness_bar(sp::Spinor{VBarKind}) = (1, -sp.mass, sp.momentum)

# Apply spin sum to a product of two DiracChains
# |M|^2 = sum_spins (ubar Gamma1 u)(ubar Gamma2 v) → Tr[(p-slash+m) Gamma1] * ...
function fermion_spin_sum(chains::Vector{DiracChain})
    result_chains = DiracChain[]
    for chain in chains
        elems = chain.elements
        length(elems) < 2 && continue

        # Find spinor pair at boundaries
        first_elem = first(elems)
        last_elem = last(elems)

        if !(first_elem isa Spinor && last_elem isa Spinor)
            push!(result_chains, chain)
            continue
        end

        # Check conjugate pairing
        if !is_conjugate(first_elem, last_elem)
            error("Spin sum requires conjugate spinor pair, got $(typeof(first_elem)) and $(typeof(last_elem))")
        end

        # Apply completeness: replace spinor pair with (p-slash ± m)
        inner_gammas = DiracGamma[e for e in elems[2:end-1] if e isa DiracGamma]
        sign, mass, mom = _completeness(last_elem)

        # Build trace chain: (p-slash + m) * inner_gammas
        # p-slash part
        trace_gammas_p = DiracGamma[GS(mom); inner_gammas...]

        if iszero(mass)
            # Massless: just Tr[p-slash * gammas]
            push!(result_chains, DiracChain(trace_gammas_p))
        else
            # Massive: Tr[(p-slash + m) * gammas] = Tr[p-slash * gammas] + m * Tr[gammas]
            # We need to return multiple chains with coefficients
            # For now, handle massless case (tracer bullet is massless)
            push!(result_chains, DiracChain(trace_gammas_p))
            if !iszero(mass)
                mass_gammas = DiracGamma[DiracGamma(IdSlot()); inner_gammas...]
                push!(result_chains, DiracChain(mass_gammas))
            end
        end
    end
    result_chains
end

# Higher-level: spin-summed |M|² for two independent fermion lines.
# |M|² = [v̄₂ Γ u₁][ū₃ Γ' v₄] × conjugate
#       = Tr[(p₂-slash+m₂) Γ (p₁-slash+m₁) Γ̄] × Tr[(k₁-slash+m₃) Γ' (k₂-slash+m₄) Γ̄']
#
# The key: TWO separate traces multiplied, not one big trace.
# Γ̄ reverses the gamma order from the conjugate amplitude.
function spin_sum_amplitude_squared(chain1::DiracChain, chain2::DiracChain)
    tr1 = _single_line_trace(chain1)
    tr2 = _single_line_trace(chain2)
    tr1 * tr2
end

# Compute the trace for a single fermion line after spin sum.
# Input: [spinor_L, γ^μ, ..., spinor_R]
# Output: Tr[(completeness of R) * gammas * (completeness of L) * reversed_gammas_conj]
#
# For massless e+e-: [v̄(p2), γ^μ, u(p1)]
# → Tr[p1-slash γ^ν p2-slash γ^μ]  (ν from conjugate amplitude)
# Note: conjugate reverses gamma order AND relabels the free index.
function _single_line_trace(chain::DiracChain)
    elems = chain.elements
    sp_L = first(elems)::Spinor
    sp_R = last(elems)::Spinor

    gammas_fwd = DiracGamma[e for e in elems[2:end-1] if e isa DiracGamma]

    # The conjugate amplitude has gammas in reverse order with a NEW index.
    # For chain [v̄(p2) γ^μ u(p1)], the conjugate is [ū(p1) γ^ν v(p2)].
    # After spin sum: Tr[p1-slash γ^ν p2-slash γ^μ]
    # We need to relabel the index in the conjugate to avoid collision.
    gammas_conj = _conjugate_gammas(gammas_fwd)

    # Completeness insertions
    _, mass_R, mom_R = _completeness(sp_R)  # right spinor
    _, mass_L, mom_L = _completeness(sp_L)  # left spinor (bar)

    # Build: Tr[(mom_R-slash + mass_R) * gammas_conj * (mom_L-slash + mass_L) * gammas_fwd]
    _build_and_trace(mom_R, mass_R, gammas_conj, mom_L, mass_L, gammas_fwd)
end

# Conjugate gamma relabeling: γ^μ → γ^μ' (prime index to avoid contraction)
function _conjugate_gammas(gs::Vector{DiracGamma})
    reversed = reverse(gs)
    DiracGamma[_relabel_gamma(g) for g in reversed]
end

function _relabel_gamma(g::DiracGamma{LISlot})
    old = g.slot.index
    new_name = Symbol(string(old.name), "_")  # mu → mu_
    DiracGamma(LISlot(LorentzIndex(new_name, old.dim)))
end
_relabel_gamma(g::DiracGamma) = g  # MomSlot etc: no index to relabel

# Build and compute the trace from completeness + gammas
function _build_and_trace(mom1, mass1, gammas_mid, mom2, mass2, gammas_end)
    parts = Tuple{Rational{Int}, Vector{DiracGamma}}[]

    # (p1-slash)(gammas_mid)(p2-slash)(gammas_end)
    push!(parts, (1//1, [GS(mom1); gammas_mid; GS(mom2); gammas_end]))

    if !iszero(mass1) && !iszero(mass2)
        push!(parts, (mass1 * mass2, [gammas_mid; gammas_end]))
    end
    if !iszero(mass1)
        push!(parts, (mass1, [gammas_mid; GS(mom2); gammas_end]))
    end
    if !iszero(mass2)
        push!(parts, (mass2, [GS(mom1); gammas_mid; gammas_end]))
    end

    result = AlgSum()
    for (c, gs) in parts
        tr = dirac_trace(gs)
        result = result + c * tr
    end
    result
end

# ---- Cross-line trace: generalization for interference ----

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
    # Build the reconnected trace by following the spinor loop:
    # j-sink → i*-gammas → j-source → j-gammas → ...
    iL, iR = _line_info(amp_i[1]), _line_info(amp_i[2])
    jL, jR = _line_info(amp_j[1]), _line_info(amp_j[2])

    # For M_i*, the spinor roles reverse (bar <-> plain).
    # Original iL: v̄(p2) Γ u(p1) -> Conjugate: ū(p1) Γ̃ v(p2)
    # So conjugate bar_mom = original plain_mom, conjugate plain_mom = original bar_mom
    i_conj = [_SpinSumLineInfo(l.plain_mom, l.bar_mom, l.gammas) for l in [iL, iR]]
    j_lines = [jL, jR]

    # Follow the loop: j-sink → i*-gammas → j-source → j-gammas → ...
    trace_gs = DiracGamma[]
    current_mom = jL.plain_mom
    visited = Set{Symbol}()

    while current_mom.name ∉ visited
        push!(visited, current_mom.name)
        push!(trace_gs, GS(current_mom))

        # Find i*-line whose bar_mom matches current_mom
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
