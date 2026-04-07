# Numerical evaluation of Passarino-Veltman scalar loop integrals.
# MS-bar scheme: UV poles (1/epsilon) subtracted. Returns finite part only.
# All functions return ComplexF64 (loop integrals are generically complex).
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf (scalar integrals A₀, B₀, C₀, D₀)
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf (PV reduction, Eqs. 4.18-4.23, App. B)
#
# B0 uses QuadGK adaptive quadrature on the Feynman parameter integral —
# correct by construction, no hand-derived analytical special cases needed.

using QuadGK: quadgk

# ---- Dispatch on PaVe{N} ----

evaluate(pv::PaVe{1}; mu2::Float64 = 1.0) = _eval_A(pv; mu2)
evaluate(pv::PaVe{2}; mu2::Float64 = 1.0) = _eval_B(pv; mu2)
evaluate(pv::PaVe{3}; mu2::Float64 = 1.0) = _eval_C(pv; mu2)
evaluate(pv::PaVe{4}; mu2::Float64 = 1.0) = _eval_D(pv; mu2)
evaluate(pv::PaVe{N}; mu2::Float64 = 1.0) where {N} =
    error("evaluate not yet implemented for PaVe{$N}")

# ---- A0: scalar 1-point function ----
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.21)
# "A₀(m) = m²(Δ - ln(m²/μ²) + 1) + O(D-4)"
# In MS-bar (Δ=0): A₀(m²) = m²(1 - ln(m²/μ²))
# A0(0) = 0  [massless tadpole vanishes in dim-reg]

function _eval_A(pv::PaVe{1}; mu2::Float64)::ComplexF64
    m2 = pv.masses[1]
    m2 == 0.0 && return 0.0 + 0.0im
    complex(m2 * (1.0 - log(m2 / mu2)))
end

# ---- C0: scalar 3-point function ----
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Eq. (5.2)
# "C = iπ² ∫₀¹ dx ∫₀^{1-x} dy [ax² + by² + cxy + dx + ey + f]^{-1}"
#
# In Minkowski convention (matching FeynCalc / our PaVe normalization):
# C₀(p₁₀, p₁₂, p₂₀, m₀², m₁², m₂²) = -∫₀¹ dx ∫₀^{1-x} dy / Δ(x,y)
# Δ(x,y) = -x(1-x-y)p₁₀ - y(1-x-y)p₂₀ - xy·p₁₂ + (1-x-y)m₀² + x·m₁² + y·m₂²
#
# UV-finite: no μ² dependence. The -iε prescription is implemented via
# complex(Δ, -ε), matching the B₀ convention.

function _eval_C(pv::PaVe{3}; mu2::Float64)::ComplexF64
    isempty(pv.indices) || return _eval_C_tensor(pv; mu2)
    _C0_analytical(pv.invariants[1], pv.invariants[2], pv.invariants[3],
                   pv.masses[1], pv.masses[2], pv.masses[3])
end

function _C0_quadgk(p10::Float64, p12::Float64, p20::Float64,
                    m02::Float64, m12::Float64, m22::Float64)::ComplexF64
    # rtol=1e-10 (vs 1e-12 for B₀): nested quadgk compounds the tolerance,
    # giving ~1e-10 overall accuracy, which is sufficient for C₀.
    val, _ = quadgk(0.0, 1.0; rtol = 1e-10) do x
        inner, _ = quadgk(0.0, 1.0 - x; rtol = 1e-10) do y
            z = 1.0 - x - y
            delta = -x * z * p10 - y * z * p20 - x * y * p12 +
                     z * m02 + x * m12 + y * m22
            -1.0 / complex(delta, -1e-30)
        end
        inner
    end
    val
end

# ---- C1, C2: tensor 3-point coefficients via Passarino-Veltman reduction ----
# Ref: refs/papers/PassarinoVeltman1979_NuclPhysB160.pdf, Sect. 4
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.6)-(4.8)
#
# C^μ = k₁^μ C₁ + k₂^μ C₂  (Lorentz decomposition of vector 3-point function)
#
# PV identity: 2l·kⱼ = Dⱼ+mⱼ²-D₀-m₀²-kⱼ²  cancels one propagator, yielding:
#   2(k₁·C) = B₀(p₂₀,m₀²,m₂²) - B₀(p₁₂,m₁²,m₂²) + (m₁²-m₀²-p₁₀)C₀
#   2(k₂·C) = B₀(p₁₀,m₀²,m₁²) - B₀(p₁₂,m₁²,m₂²) + (m₂²-m₀²-p₂₀)C₀
#
# Gram matrix: Gᵢⱼ = kᵢ·kⱼ  →  [C₁,C₂] = G⁻¹ [R₁,R₂]

function _eval_C_tensor(pv::PaVe{3}; mu2::Float64)::ComplexF64
    p10, p12, p20 = pv.invariants
    m02, m12, m22 = pv.masses
    c0 = _C0_analytical(p10, p12, p20, m02, m12, m22)
    b0_02 = _B0(p20, m02, m22; mu2)   # ∫1/(D₀D₂) → B₀(k₂², m₀², m₂²)
    b0_12 = _B0(p12, m12, m22; mu2)   # ∫1/(D₁D₂) → B₀((k₁-k₂)², m₁², m₂²)
    b0_01 = _B0(p10, m02, m12; mu2)   # ∫1/(D₀D₁) → B₀(k₁², m₀², m₁²)
    R1 = 0.5 * (b0_02 - b0_12 + (m12 - m02 - p10) * c0)
    R2 = 0.5 * (b0_01 - b0_12 + (m22 - m02 - p20) * c0)
    g11 = p10;  g12 = (p10 + p20 - p12) / 2.0;  g22 = p20
    det_G = g11 * g22 - g12^2
    abs(det_G) < 1e-30 && error("Gram matrix singular at p10=$p10, p12=$p12, p20=$p20")
    C1_val = (g22 * R1 - g12 * R2) / det_G
    C2_val = (-g12 * R1 + g11 * R2) / det_G
    pv.indices == [1] && return C1_val
    pv.indices == [2] && return C2_val
    error("PaVe{3} tensor indices $(pv.indices) not yet implemented")
end

# ---- D0: scalar 4-point function ----
# UV-finite. No μ² dependence. Dispatches to d0_collier.jl.
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Eq. (5.2)

function _eval_D(pv::PaVe{4}; mu2::Float64)::ComplexF64
    isempty(pv.indices) && return _D0_evaluate(
        pv.invariants[1], pv.invariants[2], pv.invariants[3],
        pv.invariants[4], pv.invariants[5], pv.invariants[6],
        pv.masses[1], pv.masses[2], pv.masses[3], pv.masses[4])
    _eval_D_tensor(pv; mu2)
end
