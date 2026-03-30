# D-tensor coefficients via Passarino-Veltman reduction.
#
# D^μ = k₁^μ D₁ + k₂^μ D₂ + k₃^μ D₃
#
# where k₁, k₂, k₃ are accumulated external momenta:
#   k₁² = p10, k₂² = p20, k₃² = p30
#
# PV identity: 2l·kⱼ = Dⱼ + mⱼ² - D₀ - m₀² - kⱼ²
# cancels one propagator, yielding C₀ integrals.
#
# [D₁, D₂, D₃] = G⁻¹ [R₁, R₂, R₃]
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.7)-(4.9), (4.13)-(4.18)
# Ref: refs/papers/PassarinoVeltman1979_NuclPhysB160.pdf, Sect. 4

function _eval_D_tensor(pv::PaVe{4}; mu2::Float64)::ComplexF64
    p10, p12, p23, p30, p20, p13 = pv.invariants
    m02, m12, m22, m32 = pv.masses

    # Scalar D₀
    d0 = _D0_evaluate(p10, p12, p23, p30, p20, p13, m02, m12, m22, m32)

    # C₀ integrals with each propagator removed from the box:
    # C₀^{(0)}: remove prop 0 → triangle(1,2,3), invariants (p12,p23,p13)
    c0_0 = _C0_analytical(p12, p23, p13, m12, m22, m32)
    # C₀^{(1)}: remove prop 1 → triangle(0,2,3), invariants (p20,p23,p30)
    c0_1 = _C0_analytical(p20, p23, p30, m02, m22, m32)
    # C₀^{(2)}: remove prop 2 → triangle(0,1,3), invariants (p10,p13,p30)
    c0_2 = _C0_analytical(p10, p13, p30, m02, m12, m32)
    # C₀^{(3)}: remove prop 3 → triangle(0,1,2), invariants (p10,p12,p20)
    c0_3 = _C0_analytical(p10, p12, p20, m02, m12, m22)

    # RHS: R_j = ½[C₀^{(j)} - C₀^{(0)} + f_j D₀]
    # f_j = m_j² - m₀² - k_j²  (k₁²=p10, k₂²=p20, k₃²=p30)
    R1 = 0.5 * (c0_1 - c0_0 + (m12 - m02 - p10) * d0)
    R2 = 0.5 * (c0_2 - c0_0 + (m22 - m02 - p20) * d0)
    R3 = 0.5 * (c0_3 - c0_0 + (m32 - m02 - p30) * d0)

    # Gram matrix: G_{ij} = k_i · k_j (accumulated momenta)
    g11 = p10;                         g12 = (p10 + p20 - p12) / 2.0
    g13 = (p10 + p30 - p13) / 2.0;    g22 = p20
    g23 = (p20 + p30 - p23) / 2.0;    g33 = p30

    # Solve 3×3 system via Cramer's rule
    det_G = g11 * (g22 * g33 - g23^2) -
            g12 * (g12 * g33 - g23 * g13) +
            g13 * (g12 * g23 - g22 * g13)
    abs(det_G) < 1e-30 && error("D-tensor Gram matrix singular (det=$det_G)")

    # Cofactors for each column
    D1_val = ((g22*g33 - g23^2) * R1 + (g13*g23 - g12*g33) * R2 + (g12*g23 - g13*g22) * R3) / det_G
    D2_val = ((g23*g13 - g12*g33) * R1 + (g11*g33 - g13^2) * R2 + (g12*g13 - g11*g23) * R3) / det_G
    D3_val = ((g12*g23 - g13*g22) * R1 + (g13*g12 - g11*g23) * R2 + (g11*g22 - g12^2) * R3) / det_G

    pv.indices == [1] && return D1_val
    pv.indices == [2] && return D2_val
    pv.indices == [3] && return D3_val
    error("PaVe{4} tensor indices $(pv.indices) not yet implemented")
end
