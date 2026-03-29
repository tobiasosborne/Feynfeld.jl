# SU(N) colour trace: always returns AlgSum.
# Evaluates immediately for concrete N. Default N=3 (QCD).
#
# Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eqs. (2.15)-(2.16)
# Cross-check: refs/FeynCalc/Tests/SUN/SUNTrace.test
# Key identities:
#   "Tr(1) = N"                                        [Eq. (2.15)]
#   "Tr(T^a) = 0"                                      [Eq. (2.15)]
#   "Tr(T^a T^b) = T_F δ^{ab}  where T_F = 1/2"       [Eq. (2.16)]
#   Tr(T^a T^b T^c) = (1/4)(d^{abc} + i·f^{abc})
#   n >= 4: recursive via T^a T^b = δ^{ab}/(2N) + (1/2)(d+if)T^f
#
# The recursion tracks real and imaginary parts separately:
#   trace = real_part + i × imag_part
# so that i² = -1 is handled correctly for products of f-terms.
# The public API returns the real part (sufficient for |M|²).

function colour_trace(generators::Vector{SUNT}; N::Int=3)
    re, _im = _colour_trace_ri(generators, N)
    re
end

colour_trace(chain::ColourChain; N::Int=3) = colour_trace(chain.elements; N=N)

# Returns (real::AlgSum, imag::AlgSum) where trace = real + i*imag.
function _colour_trace_ri(gs::Vector{SUNT}, N::Int)
    n = length(gs)

    # Tr(1) = N
    n == 0 && return (alg(N), AlgSum())

    # Tr(T^a) = 0
    n == 1 && return (AlgSum(), AlgSum())

    # Tr(T^a T^b) = (1/2) δ^{ab}
    if n == 2
        a, b = gs[1].adj, gs[2].adj
        return ((1//2) * alg(SUNDelta(a, b)), AlgSum())
    end

    # Tr(T^a T^b T^c) = (1/4)(d^{abc} + i·f^{abc})
    if n == 3
        a, b, c = gs[1].adj, gs[2].adj, gs[3].adj
        return ((1//4) * alg(SUND(a, b, c)), (1//4) * alg(SUNF(a, b, c)))
    end

    # n >= 4: recursive reduction
    # T^a T^b = δ^{ab}/(2N) + (1/2)(d^{abf} + i·f^{abf}) T^f
    a, b = gs[1].adj, gs[2].adj
    rest = gs[3:end]
    f_idx = _fresh_adj()

    # First term: δ^{ab}/(2N) × Tr(rest)
    re_rest, im_rest = _colour_trace_ri(rest, N)
    t1_re = (1 // (2*N)) * alg(SUNDelta(a, b)) * re_rest
    t1_im = (1 // (2*N)) * alg(SUNDelta(a, b)) * im_rest

    # Second term: (1/2)(d + if) × Tr(T^f · rest)
    # (d + if)(re_sub + i·im_sub) = (d·re - f·im) + i(d·im + f·re)
    re_sub, im_sub = _colour_trace_ri([SUNT(f_idx); rest], N)
    d = (1//2) * alg(SUND(a, b, f_idx))
    f = (1//2) * alg(SUNF(a, b, f_idx))

    t2_re = d * re_sub + (-1//1) * f * im_sub
    t2_im = d * im_sub + f * re_sub

    (t1_re + t2_re, t1_im + t2_im)
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
