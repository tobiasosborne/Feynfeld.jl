# Feynfeld.jl — SU(N) colour algebra simplification
#
# Key identities:
#   δ^{aa} = N²-1 (adjoint trace)
#   δ_F^{ii} = N (fundamental trace)
#   f^{acd} f^{bcd} = N δ^{ab}
#   d^{acd} d^{bcd} = (N²-4)/N δ^{ab}
#   T^a T^a = CF (Casimir in chain)
#   Tr(T^a T^b) = (1/2) δ^{ab}
#
# Ref: FeynCalc SUN/SUNSimplify.m

export sun_simplify, delta_trace, contract_ff, contract_dd, contract_fd
export contract_ff_full, contract_dd_full

"""
    sun_simplify(expr) -> simplified

Apply SU(N) colour algebra simplification rules.
Currently handles delta contractions and structure constant identities.
"""
function sun_simplify end

# ── Delta trace ──────────────────────────────────────────────────────

"""δ^{aa} = N²-1 (adjoint), δ_F^{ii} = N (fundamental)."""
function delta_trace(d::SUNDelta)
    d.a == d.b && return :($(SUNN)^2 - 1)
    d
end

function delta_trace(d::SUNFDelta)
    d.i == d.j && return SUNN
    d
end

# ── Structure constant contractions ──────────────────────────────────

"""
    contract_ff(f1::SUNF, f2::SUNF) -> result

Contract f^{acd} f^{bcd} = N δ^{ab}.
Returns nothing if no contraction possible.
"""
function contract_ff(f1::SUNF, f2::SUNF)
    # Find two shared indices
    shared, free1, free2 = _shared_sun_indices(
        [f1.a, f1.b, f1.c], [f2.a, f2.b, f2.c])
    shared == 2 || return nothing
    # f^{acd} f^{bcd} = N δ^{ab}
    (SUNN, SUNDelta(free1[1], free2[1]))
end

"""
    contract_dd(d1::SUND, d2::SUND) -> result

Contract d^{acd} d^{bcd} = (N²-4)/N δ^{ab}.
"""
function contract_dd(d1::SUND, d2::SUND)
    shared, free1, free2 = _shared_sun_indices(
        [d1.a, d1.b, d1.c], [d2.a, d2.b, d2.c])
    shared == 2 || return nothing
    (:(($(SUNN)^2 - 4) / $(SUNN)), SUNDelta(free1[1], free2[1]))
end

"""
    contract_fd(f::SUNF, d::SUND) -> result

f^{abc} d^{abd} = 0 (always).
"""
function contract_fd(f::SUNF, d::SUND)
    shared, _, _ = _shared_sun_indices(
        [f.a, f.b, f.c], [d.a, d.b, d.c])
    shared >= 2 && return 0
    nothing
end

# ── Full contraction identities ──────────────────────────────────────

"""f^{abc} f^{abc} = 2 N (N²-1) / (2N) = 2 CA² CF."""
function contract_ff_full(f1::SUNF, f2::SUNF)
    shared, _, _ = _shared_sun_indices(
        [f1.a, f1.b, f1.c], [f2.a, f2.b, f2.c])
    shared == 3 || return nothing
    :(2 * $(CA)^2 * $(CF))
end

"""d^{abc} d^{abc} = (N⁴-6N²+18)/(2N²) ... simplified: -2(4-CA²)CF"""
function contract_dd_full(d1::SUND, d2::SUND)
    shared, _, _ = _shared_sun_indices(
        [d1.a, d1.b, d1.c], [d2.a, d2.b, d2.c])
    shared == 3 || return nothing
    :(-2 * (4 - $(CA)^2) * $(CF))
end

# ── Helpers ──────────────────────────────────────────────────────────

"""Count shared SUNIndex names between two lists. Returns (n_shared, free1, free2)."""
function _shared_sun_indices(list1::Vector{SUNIndex}, list2::Vector{SUNIndex})
    free2 = copy(list2)
    free1 = SUNIndex[]
    shared = 0
    for a in list1
        matched = false
        for (j, b) in enumerate(free2)
            if a.name == b.name
                shared += 1
                deleteat!(free2, j)
                matched = true
                break
            end
        end
        matched || push!(free1, a)
    end
    (shared, free1, free2)
end
