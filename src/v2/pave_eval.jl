# Numerical evaluation of Passarino-Veltman scalar loop integrals.
# MS-bar scheme: UV poles (1/epsilon) subtracted. Returns finite part only.
# All functions return ComplexF64 (loop integrals are generically complex).
# Reference: 't Hooft & Veltman (1979), Denner (1993).
#
# B0 uses QuadGK adaptive quadrature on the Feynman parameter integral —
# correct by construction, no hand-derived analytical special cases needed.

using QuadGK: quadgk

# ---- Dispatch on PaVe{N} ----

evaluate(pv::PaVe{1}; mu2::Float64 = 1.0) = _eval_A(pv; mu2)
evaluate(pv::PaVe{2}; mu2::Float64 = 1.0) = _eval_B(pv; mu2)
evaluate(pv::PaVe{N}; mu2::Float64 = 1.0) where {N} =
    error("evaluate not yet implemented for PaVe{$N}")

# ---- A0: scalar 1-point function ----
# A0(m²) = m²(1 - ln(m²/μ²))  [MS-bar finite part, closed form]
# A0(0) = 0  [massless tadpole vanishes in dim-reg]

function _eval_A(pv::PaVe{1}; mu2::Float64)::ComplexF64
    m2 = pv.masses[1]
    m2 == 0.0 && return 0.0 + 0.0im
    complex(m2 * (1.0 - log(m2 / mu2)))
end

# ---- B0: scalar 2-point function ----
# B0(p², m₀², m₁²) = Δ - ∫₀¹ dx ln[f(x)/μ²]
# where f(x) = (1-x)m₀² + xm₁² - x(1-x)p² - iε
# In MS-bar (Δ=0): B0 = -∫₀¹ dx ln[f(x)/μ²]
#
# We evaluate this directly via QuadGK. The only special cases are
# where the integrand has a log singularity at x=0 or x=1 (zero masses).

function _eval_B(pv::PaVe{2}; mu2::Float64)::ComplexF64
    isempty(pv.indices) || return _eval_B_tensor(pv; mu2)
    _B0(pv.invariants[1], pv.masses[1], pv.masses[2]; mu2)
end

function _B0(p2::Float64, m02::Float64, m12::Float64; mu2::Float64)::ComplexF64
    # Both masses zero: closed form (integrand has log(x)+log(1-x) singularity)
    if m02 == 0.0 && m12 == 0.0
        return _B0_both_massless(p2; mu2)
    end
    # One mass zero: split off the log singularity analytically
    if m02 == 0.0
        return _B0_one_massless(p2, m12; mu2)
    elseif m12 == 0.0
        # B0 is symmetric in m0, m1 (Feynman parameter x → 1-x)
        return _B0_one_massless(p2, m02; mu2)
    end
    # General case (both masses nonzero): direct quadrature, no singularities
    _B0_quadgk(p2, m02, m12; mu2)
end

# Both massless: B0(p², 0, 0) = 2 - ln(-p²/μ² - iε)  [closed form]
function _B0_both_massless(p2::Float64; mu2::Float64)::ComplexF64
    p2 == 0.0 && error("B0(0, 0, 0) is IR divergent")
    # For p² > 0 (timelike): ln(-p²-iε) = ln(p²) - iπ
    2.0 - log(complex(-p2 / mu2, -1e-30))
end

# One mass zero: B0(p², 0, m²). The integrand has a log(x) singularity at x=0.
# Split: ∫₀¹ ln(x·g(x)/μ²) dx = ∫₀¹ ln(x) dx + ∫₀¹ ln(g(x)/μ²) dx = -1 + regular
# where g(x) = m² - (1-x)p²
function _B0_one_massless(p2::Float64, m2::Float64; mu2::Float64)::ComplexF64
    if p2 == 0.0
        return complex(1.0 - log(m2 / mu2))
    end
    # f(x) = xm² - x(1-x)p² = x[m² - (1-x)p²] = x · g(x)
    # B0 = -∫₀¹ ln(x·g(x)/μ²) dx = -∫₀¹ ln(x) dx - ∫₀¹ ln(g(x)/μ²) dx
    #     = 1 - ∫₀¹ ln(g(x)/μ²) dx    [since ∫₀¹ ln(x) dx = -1]
    # g(x) = m² - (1-x)p² is smooth on [0,1], no singularity → use quadgk
    integral, _ = quadgk(0.0, 1.0; rtol = 1e-12) do x
        gx = m2 - (1.0 - x) * p2
        real(log(complex(gx / mu2, -1e-30)))
    end
    complex(1.0 - integral)
end

# General case: both masses nonzero. f(x) = (1-x)m₀² + xm₁² - x(1-x)p² > 0
# for spacelike p² < 0 with positive masses, or below threshold.
# Above threshold, f(x) can become negative → log gets imaginary part.
function _B0_quadgk(p2::Float64, m02::Float64, m12::Float64; mu2::Float64)::ComplexF64
    val, _ = quadgk(0.0, 1.0; rtol = 1e-12) do x
        f = (1.0 - x) * m02 + x * m12 - x * (1.0 - x) * p2
        -log(complex(f / mu2, -1e-30))
    end
    val
end

# ---- Tensor B-functions via Passarino-Veltman reduction ----
# B1 = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)  [Denner Eq. B.7]

function _eval_B_tensor(pv::PaVe{2}; mu2::Float64)::ComplexF64
    p2, m02, m12 = pv.invariants[1], pv.masses[1], pv.masses[2]
    if pv.indices == [1]
        p2 == 0.0 && error("B1 at p²=0 requires special treatment")
        a0_m1 = evaluate(A0(m12); mu2)
        a0_m0 = evaluate(A0(m02); mu2)
        b0 = _B0(p2, m02, m12; mu2)
        return (a0_m1 - a0_m0 - (p2 + m02 - m12) * b0) / (2.0 * p2)
    end
    error("PaVe{2} tensor indices $(pv.indices) not yet implemented")
end
