# Spiral 5: QED vertex correction — anomalous magnetic moment (g-2)
#
# Ground truth: F₂(0) = α/(2π)  (Schwinger, Phys. Rev. 73 (1948) 416)
# Cross-check: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/El-GaEl.m
#              (line 210: knownResult = AlphaFS/(2*Pi))
#
# C₀ scalar 3-point function:
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Eq. (5.2)
# "C = iπ² ∫₀¹ dx ∫₀^{1-x} dy [ax² + by² + cxy + dx + ey + f]^{-1}"
#
# C₁, C₂ via Passarino-Veltman tensor reduction:
# Ref: refs/papers/PassarinoVeltman1979_NuclPhysB160.pdf, Sect. 4
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.6)-(4.8)

using Test
using Feynfeld

@testset "Spiral 5: QED vertex correction (g-2)" begin

    # ----------------------------------------------------------------
    # C₀ scalar 3-point function
    # ----------------------------------------------------------------
    @testset "C₀ evaluation" begin

        # Test 1: All momenta zero, equal masses → analytical result
        # Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Eq. (5.2)
        # Δ(x,y) = m² everywhere → C₀ = -∫₀¹ dx ∫₀^{1-x} dy / m² = -1/(2m²)
        m² = 2.0
        c0_val = evaluate(C0(0.0, 0.0, 0.0, m², m², m²))
        @test real(c0_val) ≈ -1.0 / (2.0 * m²) atol=1e-8
        @test imag(c0_val) ≈ 0.0 atol=1e-10

        # Test 2: All momenta zero, unit mass
        c0_unit = evaluate(C0(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
        @test real(c0_unit) ≈ -0.5 atol=1e-8

        # Test 3: Asymmetric masses, zero momenta
        # Δ = (1-x-y)m₀² + x m₁² + y m₂² — linear in x,y, always positive
        # Analytical: C₀(0,0,0,m₀²,m₁²,m₂²) = -∫∫ dy dx / [(1-x-y)m₀²+xm₁²+ym₂²]
        # For m₀²=1, m₁²=2, m₂²=3: compute numerically
        c0_asym = evaluate(C0(0.0, 0.0, 0.0, 1.0, 2.0, 3.0))
        @test isfinite(real(c0_asym))
        @test imag(c0_asym) ≈ 0.0 atol=1e-10
        # Cross-check: must be negative (integrand is positive)
        @test real(c0_asym) < 0.0

        # Test 4: Non-zero momenta (spacelike, below all thresholds)
        # p10=0.5, p12=0.3, p20=0.4 with heavy masses → real result
        c0_gen = evaluate(C0(0.5, 0.3, 0.4, 4.0, 4.0, 4.0))
        @test isfinite(real(c0_gen))
        @test imag(c0_gen) ≈ 0.0 atol=1e-9
        @test real(c0_gen) < 0.0

        # Test 5: g-2 kinematics (p₁²=p₂²=m², q²=0, photon mass λ)
        # C₀(m², 0, m², λ², m², m²): should be finite and negative for λ>0
        m²_e = 1.0; λ² = 0.01
        c0_g2 = evaluate(C0(m²_e, 0.0, m²_e, λ², m²_e, m²_e))
        @test isfinite(real(c0_g2))
        @test real(c0_g2) < 0.0

        # Test 6: One zero internal mass (exercises _B0_one_massless path in reduction)
        c0_zm = evaluate(C0(1.0, 0.5, 1.5, 0.0, 1.0, 2.0))
        @test isfinite(real(c0_zm))
        @test real(c0_zm) < 0.0

        # Test 7: Timelike momentum — imaginary part from iε prescription
        # p10 = 9.0 > (m0+m1)² = (0+1)² for m0²=0, m1²=1
        c0_tl = evaluate(C0(9.0, 0.5, 1.0, 0.0, 1.0, 1.0))
        @test isfinite(real(c0_tl))
        # Timelike p10 can produce nonzero imaginary part
        @test isfinite(imag(c0_tl))
    end

    # ----------------------------------------------------------------
    # C₁, C₂ tensor reduction
    # ----------------------------------------------------------------
    @testset "C₁, C₂ Passarino-Veltman reduction" begin

        # Test at generic kinematics (non-degenerate Gram matrix)
        p10 = 1.0; p12 = 0.5; p20 = 1.5
        m02 = 0.1; m12 = 1.0; m22 = 2.0

        c0_val = evaluate(C0(p10, p12, p20, m02, m12, m22))
        c1_val = evaluate(C1(p10, p12, p20, m02, m12, m22))
        c2_val = evaluate(C2(p10, p12, p20, m02, m12, m22))

        # Finiteness
        @test isfinite(real(c1_val))
        @test isfinite(real(c2_val))

        # Verify Gram matrix equation: G [C₁,C₂]ᵀ = [R₁,R₂]ᵀ
        # Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.6)-(4.8)
        b0_02 = evaluate(B0(p20, m02, m22))
        b0_12 = evaluate(B0(p12, m12, m22))
        b0_01 = evaluate(B0(p10, m02, m12))
        R1 = 0.5 * (b0_02 - b0_12 + (m12 - m02 - p10) * c0_val)
        R2 = 0.5 * (b0_01 - b0_12 + (m22 - m02 - p20) * c0_val)
        g11 = p10; g12 = (p10 + p20 - p12) / 2.0; g22 = p20
        lhs1 = g11 * c1_val + g12 * c2_val
        lhs2 = g12 * c1_val + g22 * c2_val
        @test abs(lhs1 - R1) < 1e-8
        @test abs(lhs2 - R2) < 1e-8

        # Test at second kinematics
        c1_b = evaluate(C1(2.0, 1.0, 3.0, 0.5, 1.5, 2.5))
        c2_b = evaluate(C2(2.0, 1.0, 3.0, 0.5, 1.5, 2.5))
        @test isfinite(real(c1_b))
        @test isfinite(real(c2_b))

        # Gram matrix singularity: all momenta zero → det(G)=0 → error
        @test_throws ErrorException evaluate(C1(0.0, 0.0, 0.0, 1.0, 1.0, 1.0))
    end

    # ----------------------------------------------------------------
    # F₂(0) = α/(2π)  (Schwinger result)
    # ----------------------------------------------------------------
    @testset "F₂(0) = α/(2π) — Schwinger result" begin

        α = 1.0 / 137.036
        expected = α / (2π)

        # Test 1: Direct 1D integral, massless photon limit (λ→0)
        # Ref: Schwinger, Phys. Rev. 73 (1948) 416
        # Ref: FeynCalc El-GaEl.m, line 210: "knownResult = AlphaFS/(2*Pi)"
        f2_direct = vertex_f2_zero(1.0, 0.0; alpha=α)
        @test f2_direct ≈ expected rtol=1e-10

        # Test 2: With small photon mass (IR regulator) — F₂ is IR-finite
        # Note: F₂(0,λ) differs from α/(2π) by O(λ/m × ln(m/λ)) corrections.
        # The exact Schwinger result holds only at λ=0.
        f2_lambda1 = vertex_f2_zero(1.0, 1e-4; alpha=α)
        f2_lambda2 = vertex_f2_zero(1.0, 1e-8; alpha=α)
        @test isfinite(f2_lambda1)
        @test isfinite(f2_lambda2)
        @test f2_lambda2 > f2_lambda1  # closer to α/(2π) as λ→0

        # Test 3: Convergence — λ²=1e-8 is within 0.1% of the exact result
        @test f2_lambda2 ≈ expected rtol=1e-3

        # Test 4: Different electron mass (F₂(0) is mass-independent in natural units)
        f2_heavy = vertex_f2_zero(100.0, 0.0; alpha=α)
        @test f2_heavy ≈ expected rtol=1e-10

        # Test 5: 2D integral at q²=0 agrees with 1D specialization
        f2_2d = vertex_f2(0.0, 1.0, 0.0; alpha=α)
        @test f2_2d ≈ expected rtol=1e-8

        # Test 6: α-scaling — F₂(0) is linear in α
        f2_double = vertex_f2_zero(1.0, 0.0; alpha=2α)
        @test f2_double ≈ 2 * expected rtol=1e-10
    end

    # ----------------------------------------------------------------
    # F₂(q²) at spacelike q² < 0
    # ----------------------------------------------------------------
    @testset "F₂(q²) at spacelike momentum transfer" begin

        α = 1.0 / 137.036
        m² = 1.0

        # F₂(q²) should be real for q² < 0 (spacelike) and decrease with |q²|
        f2_0 = vertex_f2(0.0, m², 0.0; alpha=α)
        f2_sp1 = vertex_f2(-1.0, m², 0.0; alpha=α)
        f2_sp2 = vertex_f2(-10.0, m², 0.0; alpha=α)

        @test f2_0 > f2_sp1 > f2_sp2 > 0.0  # monotonically decreasing
        @test f2_0 ≈ α / (2π) rtol=1e-8
    end
end
