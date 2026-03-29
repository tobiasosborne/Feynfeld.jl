# Spiral 4: QED electron self-energy at 1-loop
#
# Ground truth: refs/papers/Denner1993_FortschrPhys41.pdf, Appendix B
# Eqs. (B.6), (B.8): fermion self-energy in terms of PaVe functions
#
# QED limit (photon exchange only, Feynman gauge):
#   Σ(p) = p̸ Σ_V(p²) + m Σ_S(p²)
#   Σ_V(p²) = -(α/(4π)) [2 B₁(p², m², 0) + 1]    [Eq. B.6, QED terms only]
#   Σ_S(p²) = -(α/(4π)) [4 B₀(p², m², 0) - 2]    [Eq. B.8, QED terms only]
#
# Here m is the electron mass, the photon is massless (m₀²=0), and
# B₀, B₁ are Passarino-Veltman scalar functions evaluated in MS-bar
# (UV poles subtracted, finite parts only).
#
# B₁ formula: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (B.9)
# "B₁(p², m₀, m₁) = (m₁²-m₀²)/(2p²) (B₀(p²,m₀,m₁) - B₀(0,m₀,m₁)) - ½B₀(p²,m₀,m₁)"
#
# UV pole (P&S Eq. 10.41): Σ_div = (α/(4π))(1/ε)(4m - 2p̸)
# Ref: refs/FeynCalc/.../El-El.m, lines 125-128

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Spiral 4: QED electron self-energy (1-loop)" begin

    @testset "Self-energy via PaVe at off-shell p²" begin
        # Kinematics: electron mass m, off-shell momentum p² = 2m²
        m² = 1.0   # electron mass² (natural units)
        p² = 2.0   # off-shell: p² = 2m²
        μ² = 1.0   # renormalization scale

        # Evaluate PaVe functions (MS-bar finite parts)
        b0_val = evaluate(B0(p², 0.0, m²); mu2=μ²)
        b1_val = evaluate(B1(p², 0.0, m²); mu2=μ²)
        a0_val = evaluate(A0(m²); mu2=μ²)

        # Ref: Denner 1993, Eq. (B.6) QED limit
        # "Σ_V(p²) = -(α/(4π)) [2 B₁(p², m², 0) + 1]"
        sigma_V = -(2 * b1_val + 1)  # stripped of α/(4π) coupling

        # Ref: Denner 1993, Eq. (B.8) QED limit
        # "Σ_S(p²) = -(α/(4π)) [4 B₀(p², m², 0) - 2]"
        sigma_S = -(4 * b0_val - 2)  # stripped of α/(4π) coupling

        # Cross-check B₁ via our PV reduction formula:
        # B₁(p²,m₀²,m₁²) = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)
        b1_check = (a0_val - 0.0 - (p² + 0.0 - m²) * b0_val) / (2*p²)
        @test abs(b1_val - b1_check) < 1e-10

        # Verify B₀ and B₁ are finite (no NaN/Inf)
        @test isfinite(real(b0_val))
        @test isfinite(real(b1_val))
        @test isfinite(real(sigma_V))
        @test isfinite(real(sigma_S))

        # Check imaginary parts:
        # B₀(2, 0, 1) is above one-massless threshold (p² > m²).
        # Im(B₀) = π(p²-m²)/p² = π/2, so self-energies are complex.
        # Σ_V = -(2B₁ + 1), Σ_S = -(4B₀ - 2)
        @test imag(sigma_V) ≈ π / 4 atol=1e-8
        @test imag(sigma_S) ≈ -2π atol=1e-8

        # The self-energy evaluated at p² = m² (on-shell) determines δm and Z₂.
        # At p² = 2m² (off-shell), the finite part should be a specific value.
        # Analytical B₀(2m², 0, m²) in MS-bar:
        # B₀ = 2 - ln(m²/μ²) + (m²-p²)/p² × ln(1-p²/m²)
        # For p² = 2, m² = 1, μ² = 1: B₀ = 2 + (-1/2)ln(-1) = 2 - iπ/2
        # Note: our QuadGK evaluation computes the real part correctly;
        # the imaginary part from the iε prescription requires separate handling.
        @test real(b0_val) ≈ 2.0 atol=1e-8
    end

    @testset "Self-energy at second kinematics" begin
        # Different off-shell point: p² = 0.5 m² (spacelike-like)
        m² = 1.0; p² = 0.5; μ² = 1.0

        b0_val = evaluate(B0(p², 0.0, m²); mu2=μ²)
        b1_val = evaluate(B1(p², 0.0, m²); mu2=μ²)

        sigma_V = -(2 * b1_val + 1)
        sigma_S = -(4 * b0_val - 2)

        # For 0 < p² < m², B₀ should be real (below threshold)
        # B₀(0.5, 0, 1) = 2 - ln(1) + (1-0.5)/0.5 × ln(1-0.5) = 2 + 1 × ln(0.5)
        # = 2 + ln(0.5) = 2 - ln(2) ≈ 2 - 0.6931 ≈ 1.3069
        @test isfinite(real(b0_val))
        @test imag(b0_val) ≈ 0.0 atol=1e-10  # below threshold → real
        @test real(b0_val) ≈ 2 - log(2) atol=1e-8

        # Self-energy should be real below threshold
        @test imag(sigma_V) ≈ 0.0 atol=1e-10
        @test imag(sigma_S) ≈ 0.0 atol=1e-10
    end
end
