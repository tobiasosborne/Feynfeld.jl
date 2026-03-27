# SU(N) colour trace: always returns AlgSum.
# Evaluates immediately for concrete N. Default N=3 (QCD).
#
# Key identities:
#   Tr(1) = N
#   Tr(T^a) = 0
#   Tr(T^a T^b) = T_F δ^{ab}  where T_F = 1/2
#   Tr(T^a T^b T^c) = (1/4)(d^{abc} + i·f^{abc})
#     → real part: (1/4)d^{abc}, imaginary part: (1/4)f^{abc}
#   n ≥ 4: recursive reduction via T^a T^b = δ^{ab}/(2N) + (1/2)(d+if)T^f

function colour_trace(generators::Vector{SUNT}; N::Int=3)
    n = length(generators)
    _colour_trace_impl(generators, N)
end

colour_trace(chain::ColourChain; N::Int=3) = colour_trace(chain.elements; N=N)

function _colour_trace_impl(gs::Vector{SUNT}, N::Int)
    n = length(gs)

    # Tr(1) = N
    n == 0 && return alg(N)

    # Tr(T^a) = 0
    n == 1 && return AlgSum()

    # Tr(T^a T^b) = (1/2) δ^{ab}
    if n == 2
        a, b = gs[1].adj, gs[2].adj
        return alg(SUNDelta(a, b)) * alg(1//2)
    end

    # Tr(T^a T^b T^c) = (1/4)(d^{abc} + i·f^{abc})
    # We track only the REAL part for now (sufficient for |M|² calculations).
    # The f^{abc} term contributes to imaginary part and cancels in |M|².
    if n == 3
        a, b, c = gs[1].adj, gs[2].adj, gs[3].adj
        d_term = SUND(a, b, c)
        f_term = SUNF(a, b, c)
        # Return: (1/4) d^{abc} as the symmetric part
        # and (1/4) f^{abc} as the antisymmetric part (with i)
        # For real amplitudes, keep both — they appear in different interference terms
        return (1//4) * alg(d_term) + (1//4) * alg(f_term)
    end

    # n ≥ 4: recursive reduction
    # T^a T^b = δ^{ab}/(2N) · I + (1/2)(d^{abf} + i·f^{abf}) T^f
    # Tr(T^a T^b ... T^n) = δ^{ab}/(2N) · Tr(T^c...T^n)
    #                      + (1/2) Σ_f (d^{abf} + i·f^{abf}) Tr(T^f T^c...T^n)
    a, b = gs[1].adj, gs[2].adj
    rest = gs[3:end]
    f_idx = _fresh_adj()

    # First term: δ^{ab}/(2N) · Tr(rest)
    sub1 = _colour_trace_impl(rest, N)
    term1 = (1 // (2*N)) * alg(SUNDelta(a, b)) * sub1

    # Second term: (1/2)(d^{abf} + i·f^{abf}) · Tr(T^f · rest)
    new_chain = [SUNT(f_idx); rest]
    sub2 = _colour_trace_impl(new_chain, N)
    d_part = (1//2) * alg(SUND(a, b, f_idx)) * sub2
    f_part = (1//2) * alg(SUNF(a, b, f_idx)) * sub2

    term1 + d_part + f_part
end

# ---- Colour delta contraction ----
# δ^{aa} = N² - 1 (adjoint dimension)
# δ_{ii} = N (fundamental dimension)
function colour_delta_trace(d::SUNDelta; N::Int=3)
    d.a == d.b ? Rational{Int}(N^2 - 1) : nothing
end

function colour_delta_trace(d::FundDelta; N::Int=3)
    d.i == d.j ? Rational{Int}(N) : nothing
end

# ---- Casimirs (computed for concrete N) ----
casimir_fundamental(N::Int) = (N^2 - 1) // (2*N)  # CF = (N²-1)/(2N)
casimir_adjoint(N::Int) = Rational{Int}(N)         # CA = N
trace_normalization(N::Int) = 1//2                  # T_F = 1/2
