# Dirac trace: always returns AlgSum. No mixed return types.
# MomentumSum gammas are expanded before tracing (linearity of trace).
# gamma5 traces produce Eps (Levi-Civita) tensors.
#
# Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eqs. (2.15)-(2.18)
# "Tr[γ5 γ^a γ^b γ^c γ^d] = -4i ε^{abcd}"  [FeynCalc convention, $LeviCivitaSign=-1]
# Cross-check: refs/FeynCalc/Tests/Dirac/DiracTrace.test, IDs 11-28

function dirac_trace(chain::DiracChain)
    gs = gammas(chain)
    _expand_and_trace(gs)
end

function dirac_trace(gammas_list::Vector{<:DiracGamma})
    _expand_and_trace(gammas_list)
end

# Pre-processing pipeline: expand MomSumSlot → expand projectors → dispatch to trace
function _expand_and_trace(gs::Vector{<:DiracGamma})
    # 1. Expand MomSumSlot (linearity of trace)
    idx = findfirst(g -> g isa DiracGamma{MomSumSlot}, gs)
    if idx !== nothing
        g = gs[idx]::DiracGamma{MomSumSlot}
        result = AlgSum()
        for (c, m) in g.slot.mom.terms
            expanded = copy(gs)
            expanded[idx] = GS(m)
            result = result + c * _expand_and_trace(expanded)
        end
        return result
    end

    # 2. Expand chiral projectors: GA6 → (I + γ5)/2, GA7 → (I - γ5)/2
    pidx = findfirst(g -> g isa DiracGamma{ProjPSlot} || g isa DiracGamma{ProjMSlot}, gs)
    if pidx !== nothing
        g = gs[pidx]
        sign_g5 = g isa DiracGamma{ProjPSlot} ? 1 : -1
        # Replace with identity term + γ5 term
        gs_id = DiracGamma[gs[i] for i in eachindex(gs) if i != pidx]
        gs_g5 = DiracGamma[i == pidx ? GA5() : gs[i] for i in eachindex(gs)]
        return (1 // 2) * _expand_and_trace(gs_id) + (sign_g5 // 2) * _expand_and_trace(gs_g5)
    end

    # 3. Dispatch: γ5 present or not?
    has_g5 = any(g -> g isa DiracGamma{Gamma5Slot}, gs)
    has_g5 ? _trace_with_g5(gs) : _trace_no_g5(gs)
end

# ---- Core trace: no γ5 ----
function _trace_no_g5(gs::Vector{<:DiracGamma})
    n = length(gs)
    n == 0 && return alg(4)
    isodd(n) && return AlgSum()

    if n == 2
        p = gamma_pair(gs[1], gs[2])
        return 4 * p
    end

    # Recursive: Tr[γ^a1 ... γ^an] = Σ_{k=2}^{n} (-1)^k g^{a1,ak} Tr[rest]
    result = AlgSum()
    g1 = gs[1]
    for k in 2:n
        sign = iseven(k) ? 1 : -1
        p = gamma_pair(g1, gs[k])
        iszero(p) && continue
        rest = [gs[i] for i in 2:n if i != k]
        sub_trace = _trace_no_g5(rest)
        result = result + sign * (p * sub_trace)
    end
    result
end

# ---- Trace with γ5: anticommute γ5 to end, then apply formula ----
# γ5 anticommutes with all γ^μ: γ5 γ^μ = -γ^μ γ5
# So Tr[... γ5 ... γ^μ ...] = (-1)^k × Tr[... γ^μ ... γ5] where k = moves
function _trace_with_g5(gs::Vector{<:DiracGamma})
    # Strip γ5 and collect sign from anticommutation
    ordinary = DiracGamma[]
    g5_count = 0
    g5_sign = 1
    for (i, g) in enumerate(gs)
        if g isa DiracGamma{Gamma5Slot}
            g5_count += 1
            # Sign from anticommuting γ5 past the remaining gammas to the right
            remaining_after = length(gs) - i - count(g2 -> g2 isa DiracGamma{Gamma5Slot}, gs[i+1:end])
            g5_sign *= iseven(remaining_after) ? 1 : -1
        else
            push!(ordinary, g)
        end
    end

    # γ5² = I, so even number of γ5 → trace without γ5, odd → trace with one γ5
    if iseven(g5_count)
        return g5_sign * _trace_no_g5(ordinary)
    end

    # Odd γ5: need Tr[γ^{a1}...γ^{an} γ5]
    n = length(ordinary)

    # Tr[γ5] = 0, Tr[γ^a γ5] = 0, Tr[γ^a γ^b γ5] = 0, Tr[γ^a γ^b γ^c γ5] = 0
    n < 4 && return AlgSum()
    # Odd number of ordinary gammas + γ5 → total even but no pairing → 0
    isodd(n) && return AlgSum()

    # Base case: Tr[γ^a γ^b γ^c γ^d γ5] = -4i ε^{abcd}
    # Ref: MBD Eq. (2.18), FeynCalc convention $LeviCivitaSign = -1
    if n == 4
        eps = Eps(_slot_to_pairarg(ordinary[1]), _slot_to_pairarg(ordinary[2]),
                  _slot_to_pairarg(ordinary[3]), _slot_to_pairarg(ordinary[4]))
        # -4i × sign from anticommutation
        return g5_sign * (-4) * alg(eps)  # -4i, but i is implicit (Eps IS the imaginary structure)
    end

    # Recursive: Tr[γ^{a1}...γ^{a2n} γ5]
    # Same recursion as _trace_no_g5 but sub-traces keep γ5
    result = AlgSum()
    g1 = ordinary[1]
    for k in 2:n
        sign = iseven(k) ? 1 : -1
        p = gamma_pair(g1, ordinary[k])
        iszero(p) && continue
        rest = [ordinary[i] for i in 2:n if i != k]
        # Rest still has γ5 (append it for recursion)
        sub_trace = _trace_with_g5_pure(rest)
        result = result + sign * (p * sub_trace)
    end
    g5_sign * result
end

# Pure γ5 trace: ordinary gammas only, one implicit γ5 at the end
function _trace_with_g5_pure(gs::Vector{<:DiracGamma})
    n = length(gs)
    n < 4 && return AlgSum()
    isodd(n) && return AlgSum()

    if n == 4
        eps = Eps(_slot_to_pairarg(gs[1]), _slot_to_pairarg(gs[2]),
                  _slot_to_pairarg(gs[3]), _slot_to_pairarg(gs[4]))
        return (-4) * alg(eps)
    end

    result = AlgSum()
    g1 = gs[1]
    for k in 2:n
        sign = iseven(k) ? 1 : -1
        p = gamma_pair(g1, gs[k])
        iszero(p) && continue
        rest = [gs[i] for i in 2:n if i != k]
        sub_trace = _trace_with_g5_pure(rest)
        result = result + sign * (p * sub_trace)
    end
    result
end

# Extract PairArg from a DiracGamma slot (for building Eps)
_slot_to_pairarg(g::DiracGamma{LISlot}) = g.slot.index
_slot_to_pairarg(g::DiracGamma{MomSlot}) = g.slot.mom
