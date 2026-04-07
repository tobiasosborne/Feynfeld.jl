# B₀ scalar 2-point function: numerical evaluation.
# MS-bar scheme (Δ=0): B₀ = -∫₀¹ dx ln[f(x)/μ²]
# where f(x) = (1-x)m₀² + xm₁² - x(1-x)p² - iε
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.23)
# "B₀(p², m₀, m₁) = Δ - ∫₀¹ dx ln[(p²x²-x(p²-m₀²+m₁²)+m₀²-iε)/μ²]"
#
# Uses QuadGK adaptive quadrature on the Feynman parameter integral —
# correct by construction, no hand-derived analytical special cases needed.

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
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.23) with m₀=m₁=0
# Setting m₀=m₁=0, MS-bar (Δ=0): B₀(p²,0,0) = 2 - ln(-p²/μ²)
function _B0_both_massless(p2::Float64; mu2::Float64)::ComplexF64
    if p2 == 0.0
        # B₀(0,0,0) is IR divergent in D=4. In dim-reg (D=4-2ε), the
        # 1/ε pole is absorbed into the UV subtraction; COLLIER returns
        # the finite remainder (= 0 in MS-bar for this degenerate case).
        # Ref: Denner1993, Eq. (4.23) with all arguments zero
        return _B0_collier(0.0, 0.0, 0.0)
    end
    # For p² > 0 (timelike): ln(-p²-iε) = ln(p²) - iπ
    2.0 - log(complex(-p2 / mu2, -1e-30))
end

# COLLIER B₀ via b0_coli_: handles IR-divergent cases in dimensional regularization.
# Returns the finite part in MS-bar (UV poles subtracted).
function _B0_collier(p2::Float64, m02::Float64, m12::Float64)::ComplexF64
    _ensure_collier_init()
    ccall((:b0_coli_, _COLLIER_LIB), ComplexF64,
          (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
          ComplexF64(p2), ComplexF64(m02), ComplexF64(m12))
end

# One mass zero: B0(p², 0, m²). The integrand has a log(x) singularity at x=0.
# Split: ∫₀¹ ln(x·g(x)/μ²) dx = ∫₀¹ ln(x) dx + ∫₀¹ ln(g(x)/μ²) dx = -1 + regular
# where g(x) = m² - (1-x)p²
function _B0_one_massless(p2::Float64, m2::Float64; mu2::Float64)::ComplexF64
    if p2 == 0.0
        return complex(1.0 - log(m2 / mu2))
    end
    # Real part via complex-log quadrature (handles g(x)=0 smoothly)
    integral, _ = quadgk(0.0, 1.0; rtol = 1e-12) do x
        gx = m2 - (1.0 - x) * p2
        real(log(complex(gx / mu2, -1e-30)))
    end
    real_part = 1.0 - integral

    # Imaginary part via Kallen: λ(p², 0, m²) = (p² - m²)²
    # Im(B₀) = π(p²-m²)/p² when p² > m²
    # Ref: derived from -iε prescription, f(x)<0 region width = (p²-m²)/p²
    imag_part = p2 > m2 ? π * (p2 - m2) / p2 : 0.0

    complex(real_part, imag_part)
end

# General case: both masses nonzero.
# Below/at threshold: f(x) >= 0, original complex method (numerically stable).
# Above threshold: real part via |f| quadrature, imaginary part via Kallen function.
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.23)
function _B0_quadgk(p2::Float64, m02::Float64, m12::Float64; mu2::Float64)::ComplexF64
    # Kallen function determines threshold
    # λ = (p² - (m₀+m₁)²)(p² - (m₀-m₁)²)
    lambda = p2^2 - 2.0 * p2 * (m02 + m12) + (m02 - m12)^2

    if lambda > 0.0 && p2 > 0.0
        # Above threshold: f(x) changes sign at two roots.
        # Use |f| for real part (avoids complex log), analytical imaginary part.
        real_part, _ = quadgk(0.0, 1.0; rtol = 1e-10) do x
            f = (1.0 - x) * m02 + x * m12 - x * (1.0 - x) * p2
            -log(max(abs(f / mu2), 1e-300))
        end
        imag_part = π * sqrt(lambda) / p2
        complex(real_part, imag_part)
    else
        # Below/at threshold: f(x) >= 0, no sign change. Original method is stable.
        val, _ = quadgk(0.0, 1.0; rtol = 1e-12) do x
            f = (1.0 - x) * m02 + x * m12 - x * (1.0 - x) * p2
            -log(complex(f / mu2, -1e-30))
        end
        val
    end
end

# ---- Tensor B-functions via Passarino-Veltman reduction ----
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (B.9)
# "B₁(p², m₀, m₁) = (m₁²-m₀²)/(2p²)(B₀(p²,m₀,m₁) - B₀(0,m₀,m₁)) - ½B₀(p²,m₀,m₁)"
# Equivalent form via Eq. (4.18) PV reduction + Eq. (4.21) A₀:
# B₁ = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)

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
