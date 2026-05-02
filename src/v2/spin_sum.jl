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
            push!(result_chains, DiracChain(trace_gammas_p))
            mass_gammas = DiracGamma[DiracGamma(IdSlot()); inner_gammas...]
            push!(result_chains, DiracChain(mass_gammas))
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

# DiracExpr version: handles chiral vertices (sum of chains with coefficients).
function spin_sum_amplitude_squared(de1::DiracExpr, de2::DiracExpr)
    _single_line_trace(de1) * _single_line_trace(de2)
end

# Spin-summed trace for a DiracExpr fermion line.
# |Σ_i c_i ū Γ_i v|² = Σ_{i,j} c_i c_j × Tr[compl_R × Γ̄_j × compl_L × Γ_i]
function _single_line_trace(de::DiracExpr)
    first_chain = de.terms[1][2]
    sp_L = first(first_chain.elements)::Spinor
    sp_R = last(first_chain.elements)::Spinor
    _, mass_R, mom_R = _completeness(sp_R)
    _, mass_L, mom_L = _completeness(sp_L)

    result = AlgSum()
    for (ci, chain_i) in de.terms
        gammas_i = DiracGamma[e for e in chain_i.elements[2:end-1]]
        for (cj, chain_j) in de.terms
            gammas_j = DiracGamma[e for e in chain_j.elements[2:end-1]]
            sign_j = _conj_gamma5_sign(gammas_j)
            gammas_conj_j = _conjugate_gammas(gammas_j)
            tr = _build_and_trace(mom_R, mass_R, gammas_conj_j, mom_L, mass_L, gammas_i)
            mul_acc!(result, ci * cj, tr, sign_j)
        end
    end
    result
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
# Dirac conjugation: γ5 → -γ5, GA6 ↔ GA7 (γ^0 γ5 γ^0 = -γ5)
_relabel_gamma(::DiracGamma{ProjPSlot}) = GA7()   # (1+γ5)/2 → (1-γ5)/2
_relabel_gamma(::DiracGamma{ProjMSlot}) = GA6()   # (1-γ5)/2 → (1+γ5)/2
_relabel_gamma(g::DiracGamma) = g  # MomSlot etc: no index to relabel

# Count γ5 sign flip in conjugate: each bare γ5 contributes (-1)
function _conj_gamma5_sign(gs::Vector{DiracGamma})
    n = count(g -> g isa DiracGamma{Gamma5Slot}, gs)
    iseven(n) ? 1//1 : -1//1
end

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
        add!(result, dirac_trace(gs), c)
    end
    result
end
