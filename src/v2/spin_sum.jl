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
