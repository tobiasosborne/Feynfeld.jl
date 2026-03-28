# Spiral 7: Tree-level e⁺e⁻ → W⁺W⁻ cross-section
#
# Ground truth:
# Ref: refs/FeynCalc/FeynCalc/Examples/EW/Tree/Mathematica/AnelEl-WW.m
#      lines 193-198 (Grozin formula), line 176 (numerical plot at LEP2 energies)
# Ref: Grozin, "Using REDUCE in High Energy Physics", Chapter 5.4
# Ref: Altarelli et al., "Physics at LEP2: Vol. 1", p. 93, Fig. 2
#      (experimental comparison: σ ≈ 15-17 pb at √s = 200 GeV)
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.1)-(11.2)
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eqs. (10.22a), (10.63)
#      "sin²θ_W = 0.22348", "M_W = 80.360 GeV"

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Spiral 7: e⁺e⁻ → W⁺W⁻" begin

    # ----------------------------------------------------------------
    # EW parameters
    # ----------------------------------------------------------------
    @testset "EW parameters" begin
        @test EW_M_W ≈ 80.360 atol=0.001
        @test EW_M_Z ≈ 91.188 atol=0.001
        @test EW_SIN2_W ≈ 0.22348 atol=1e-5
        @test EW_COS2_W ≈ 1.0 - EW_SIN2_W rtol=1e-12
        # Tree-level mass relation: M_Z ≈ M_W / cos(θ_W)
        mz_tree = EW_M_W / EW_COS_W
        @test mz_tree ≈ EW_M_Z rtol=0.005  # 0.5% (radiative corrections shift this)
        # Z coupling: g_V^e = -1/2 + 2sin²θ_W
        @test EW_GV_E ≈ -0.5 + 2 * EW_SIN2_W rtol=1e-10
        @test EW_GA_E ≈ -0.5 rtol=1e-10
    end

    # ----------------------------------------------------------------
    # Massive polarization sum
    # ----------------------------------------------------------------
    @testset "Massive polarization sum" begin
        mu = LorentzIndex(:mu)
        nu = LorentzIndex(:nu)
        k = Momentum(:k)
        M2 = 1 // 1  # M² = 1

        pol = polarization_sum_massive(mu, nu, k, M2)
        # Should be -g^{μν} + k^μ k^ν / M²
        @test pol isa AlgSum
        # Contains metric tensor term and momentum term
        @test length(pol.terms) == 2
    end

    # ----------------------------------------------------------------
    # Threshold behavior
    # ----------------------------------------------------------------
    @testset "Threshold" begin
        M_W = EW_M_W
        # Below threshold: σ = 0
        @test sigma_ee_ww(160.0^2) == 0.0  # √s = 160 < 2M_W ≈ 161

        # Just above threshold: σ > 0 and small (β³ behavior)
        s_above = (2 * M_W + 1.0)^2  # √s ≈ 162 GeV
        @test sigma_ee_ww(s_above) > 0.0
        @test sigma_ee_ww(s_above) < 8.0  # small near threshold

        # Threshold is at 4M_W² in physical units
        s_thresh = 4.0 * M_W^2
        @test sigma_ee_ww(s_thresh - 1.0) == 0.0
        @test sigma_ee_ww(s_thresh + 100.0) > 0.0
    end

    # ----------------------------------------------------------------
    # Cross-section at LEP2 energies
    # ----------------------------------------------------------------
    @testset "σ at LEP2 energies" begin
        # Ref: FeynCalc AnelEl-WW.m, line 176 (plot from √s=162 to 205 GeV)
        # Ref: Altarelli:1996gh, p. 93, Fig. 2 (measured ≈ 15-17 pb at 200 GeV)
        # Tree-level is slightly higher than measured (NLO corrections ~ -10%)

        # σ(√s = 170 GeV): near threshold, rapidly rising
        sig_170 = sigma_ee_ww(170.0^2)
        @test 10.0 < sig_170 < 20.0

        # σ(√s = 190 GeV): near peak
        sig_190 = sigma_ee_ww(190.0^2)
        @test 15.0 < sig_190 < 22.0

        # σ(√s = 200 GeV): established LEP2 energy
        sig_200 = sigma_ee_ww(200.0^2)
        @test 15.0 < sig_200 < 22.0

        # Rising from threshold to ~190 GeV
        @test sigma_ee_ww(170.0^2) < sigma_ee_ww(190.0^2)

        # Slowly decreasing at high energy (unitarity: σ ~ 1/s)
        sig_500 = sigma_ee_ww(500.0^2)
        @test sig_500 < sig_200
        @test sig_500 > 0.0
    end

    # ----------------------------------------------------------------
    # High-energy behavior (gauge cancellation)
    # ----------------------------------------------------------------
    @testset "High-energy behavior" begin
        # Without triple gauge couplings, σ ~ s/M_W⁴ (grows, violates unitarity).
        # WITH gauge cancellations (SM), σ ~ 1/s at high energy.
        # Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.16)-(11.17)

        sig_500 = sigma_ee_ww(500.0^2)
        sig_1000 = sigma_ee_ww(1000.0^2)

        # σ decreases with energy: σ(1000) < σ(500)
        @test sig_1000 < sig_500

        # Scaling: σ ~ 1/s at very high energy (up to log corrections)
        # σ(1000)/σ(500) ≈ (500/1000)² = 0.25 (approximately)
        ratio = sig_1000 / sig_500
        @test 0.1 < ratio < 0.5  # broad tolerance for log corrections
    end

    # ----------------------------------------------------------------
    # Consistency checks
    # ----------------------------------------------------------------
    @testset "Consistency" begin
        # Cross-section is positive at all accessible energies
        for sqrts in [165.0, 175.0, 200.0, 300.0, 500.0, 1000.0]
            @test sigma_ee_ww(sqrts^2) > 0.0
        end

        # α-dependence: σ ∝ α² (tree-level)
        sig_std = sigma_ee_ww(200.0^2; alpha=EW_ALPHA)
        sig_2a = sigma_ee_ww(200.0^2; alpha=2*EW_ALPHA)
        # Not exactly 4× because sin²θ_W also depends on α in the full theory,
        # but in our fixed-parameter approach, σ scales with α²
        @test sig_2a / sig_std ≈ 4.0 rtol=0.01

        # Different M_W: lighter W → lower threshold, different σ
        sig_light = sigma_ee_ww(200.0^2; M_W=70.0)
        @test sig_light != sig_std
        @test sig_light > 0.0
    end

    # ----------------------------------------------------------------
    # Comparison with QED tree-level (sanity check)
    # ----------------------------------------------------------------
    @testset "Comparison with QED e⁺e⁻→μ⁺μ⁻" begin
        s = 200.0^2
        GEV2_TO_PB = 3.894e8

        # QED tree: σ = 4πα²/(3s) ≈ 0.092 pb at √s = 200 GeV
        sigma_qed = 4π * EW_ALPHA^2 / (3 * s) * GEV2_TO_PB
        sigma_ww = sigma_ee_ww(s)

        # W pair production is larger than μ pair (W coupling ~ 1/sin²θ_W)
        @test sigma_ww / sigma_qed > 5.0
        @test sigma_ww / sigma_qed < 50.0
    end
end
