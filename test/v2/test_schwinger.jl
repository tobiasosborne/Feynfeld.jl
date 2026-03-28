# Integration test: Schwinger correction for e+e- → μ+μ-.
# Validates Layer 5 (PaVe evaluation) integrated with Layer 6 (cross-section).
# Ground truth: analytical formula from P&S Ch. 6-7.

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

const α = 1 / 137.036

@testset "Schwinger correction" begin
    @testset "Analytical formula" begin
        expected = 3 * α / (4 * π) * (π^2 / 3 - 1 / 2)
        delta = schwinger_correction(; alpha = α)

        @test delta ≈ expected rtol = 1e-12
        # Sanity: O(α) correction is small but positive
        @test 0.001 < delta < 0.01
        # Known approximate value
        @test delta ≈ 0.004860 atol = 0.0001
    end

    @testset "NLO cross-section" begin
        s = 100.0^2  # (100 GeV)² — high energy
        sigma_tree = 4 * π * α^2 / (3 * s)
        sigma_nlo = sigma_nlo_ee_mumu(s; alpha = α)

        # NLO should be slightly larger than tree (positive correction)
        @test sigma_nlo > sigma_tree
        @test sigma_nlo / sigma_tree ≈ 1 + schwinger_correction(; alpha = α) rtol = 1e-10
    end

    @testset "Vacuum polarization via PaVe" begin
        # Test that our B₀-based vacuum polarization gives sensible results.
        # For a single electron loop at s >> m_e²:
        m_e2 = 0.000511^2  # (electron mass in GeV)²

        # Below threshold (s < 4m²): Π̂ should be real and small
        s_low = 0.5 * m_e2
        vp_low = vacuum_polarization(s_low, m_e2; alpha = α)
        @test abs(vp_low) < 0.01  # small below threshold

        # Well above threshold: Π̂ is negative (screening)
        # At s = 1 GeV², the running is Π̂ ≈ -(α/3π) ln(s/m_e²) ≈ -0.033
        s_high = 1.0  # 1 GeV²
        vp_high = vacuum_polarization(s_high, m_e2; alpha = α)
        @test vp_high < 0  # screening
        # Rough estimate: -(α/3π) × ln(1/m_e²) ≈ -(1/137)/(3π) × 15 ≈ -0.012
        @test -0.05 < vp_high < 0.0
    end

    @testset "Vacuum polarization known limits" begin
        # For s = 0: Π̂(0) = 0 by definition (on-shell renormalization)
        m2 = 1.0
        @test vacuum_polarization(0.0, m2; alpha = α) ≈ 0.0 atol = 1e-10

        # For s >> m² (leading log approximation):
        # Π̂(s) ≈ -(α/(3π)) × [ln(s/m²) - 5/3]
        s = 1e8  # s/m² = 10⁸ >> 1
        vp = vacuum_polarization(s, m2; alpha = α)
        leading_log = -α / (3 * π) * (log(s / m2) - 5 / 3)
        @test vp ≈ leading_log rtol = 0.05  # 5% accuracy for leading-log
    end

    @testset "PaVe evaluation spot checks at Schwinger kinematics" begin
        # B₀(s, m², m²) for vacuum polarization
        # At s = (91.2 GeV)², m_e = 0.000511 GeV:
        s = 91.2^2
        m_e2 = 0.000511^2
        b0_vp = evaluate(B0(s, m_e2, m_e2); mu2 = s)
        # Should be complex (above threshold 4m_e²)
        @test imag(b0_vp) != 0.0

        # B₀(0, 0, m²) for self-energy
        b0_se = evaluate(B0(0.0, 0.0, m_e2); mu2 = m_e2)
        @test real(b0_se) ≈ 1.0 atol = 1e-10  # 1 - ln(m²/μ²) = 1 at μ² = m²

        # A₀(m²) at μ² = m²
        @test real(evaluate(A0(m_e2); mu2 = m_e2)) ≈ m_e2 atol = 1e-20
    end
end
