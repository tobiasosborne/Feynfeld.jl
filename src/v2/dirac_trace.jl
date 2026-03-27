# Dirac trace: always returns AlgSum. No mixed return types.
# MomentumSum gammas are expanded before tracing (linearity of trace).

function dirac_trace(chain::DiracChain)
    gs = gammas(chain)
    _expand_and_trace(gs)
end

function dirac_trace(gammas_list::Vector{<:DiracGamma})
    _expand_and_trace(gammas_list)
end

# Expand any MomSumSlot gammas, then trace each expanded chain.
# Tr[... γ·(a+b) ...] = Tr[... γ·a ...] + Tr[... γ·b ...]
function _expand_and_trace(gs::Vector{<:DiracGamma})
    # Find first MomSumSlot
    idx = findfirst(g -> g isa DiracGamma{MomSumSlot}, gs)
    if idx === nothing
        return _trace_no_g5(gs)  # no MomSumSlot: trace directly
    end

    # Expand at position idx
    g = gs[idx]::DiracGamma{MomSumSlot}
    result = AlgSum()
    for (c, m) in g.slot.mom.terms
        expanded = copy(gs)
        expanded[idx] = GS(m)  # replace MomSumSlot with single Momentum
        result = result + c * _expand_and_trace(expanded)
    end
    result
end

# Core trace: no MomSumSlot gammas, no γ5.
function _trace_no_g5(gs::Vector{<:DiracGamma})
    n = length(gs)
    n == 0 && return alg(4)
    isodd(n) && return AlgSum()

    # Tr[γ^a γ^b] = 4 g^{ab}
    if n == 2
        p = gamma_pair(gs[1], gs[2])
        p === nothing && return AlgSum()
        p isa Number && return alg(4 * p)
        return 4 * alg(p)
    end

    # Recursive: Tr[γ^a1 ... γ^an] = Σ_{k=2}^{n} (-1)^k g^{a1,ak} Tr[rest]
    result = AlgSum()
    g1 = gs[1]
    for k in 2:n
        sign = iseven(k) ? 1 : -1
        p = gamma_pair(g1, gs[k])
        p === nothing && continue

        rest = [gs[i] for i in 2:n if i != k]
        sub_trace = _trace_no_g5(rest)

        if p isa Number
            result = result + (sign * p) * sub_trace
        else
            result = result + sign * (alg(p) * sub_trace)
        end
    end
    result
end
