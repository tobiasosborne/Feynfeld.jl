# Spiral 6: Running fine-structure constant α(q²) from vacuum polarization.
#
# Ground truth:
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 3 (renormalization)
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eq. (10.11)-(10.12)
#      "α⁻¹(0) = 137.035999178(8)"
#      "Δα_had^(5)(M_Z) = 0.02783 ± 0.00006"
#      Δα_lep(M_Z²) ≈ 0.03150  (perturbative, 3-flavor lepton loops)
#
# Vacuum polarization formula (for one fermion, charge Q, mass m, N_c colors):
#   Π̂_f(q²) = -(2α Q² N_c / π) ∫₀¹ dx x(1-x) ln[1 - x(1-x)q²/m²]
#
# Running coupling:
#   α(q²) = α / (1 - Δα(q²)),  Δα = -Σ_f Π̂_f

using Test
using Feynfeld

const α = 1.0 / 137.036
const M_Z = 91.1876  # Z boson mass in GeV (PDG)

@testset "Spiral 6: Running α(q²)" begin

    # ----------------------------------------------------------------
    # SM fermion table
    # ----------------------------------------------------------------
    @testset "SM fermion table" begin
        @test length(SM_LEPTONS) == 3
        @test length(SM_QUARKS) == 6
        @test length(SM_FERMIONS) == 9

        # Electron is lightest
        @test SM_LEPTONS[1].mass2 < SM_LEPTONS[2].mass2 < SM_LEPTONS[3].mass2
        # Color factors
        @test all(f -> f.nc == 1, SM_LEPTONS)
        @test all(f -> f.nc == 3, SM_QUARKS)
    end

    # ----------------------------------------------------------------
    # Δα: vacuum polarization sum
    # ----------------------------------------------------------------
    @testset "Δα leptonic contribution" begin
        # Leptonic Δα at M_Z² is known precisely: 0.03150
        # Ref: PDG "Electroweak Model and Constraints on New Physics"
        da_lep = real(delta_alpha(M_Z^2; alpha=α, fermions=SM_LEPTONS))
        @test da_lep ≈ 0.0315 rtol=0.02  # 2% tolerance

        # At q² = 0: Δα = 0 (on-shell renormalization)
        @test real(delta_alpha(0.0; alpha=α)) ≈ 0.0 atol=1e-10

        # Δα should be positive (screening → α increases with energy)
        @test da_lep > 0.0

        # Electron-only contribution at M_Z
        da_e = real(delta_alpha(M_Z^2; alpha=α, fermions=SM_LEPTONS[1:1]))
        @test 0.01 < da_e < 0.02  # electron dominates leptonic running
    end

    @testset "Δα total (perturbative)" begin
        # Total Δα at M_Z² (perturbative quarks + leptons)
        # PDG value: Δα(M_Z²) ≈ 0.0590 (combined lep + had)
        # Perturbative quarks differ from data-driven Δα_had by ~10-20%
        da_total = real(delta_alpha(M_Z^2; alpha=α))
        @test 0.04 < da_total < 0.08  # broad tolerance for perturbative quarks

        # Hadronic > leptonic (more flavors, color factor 3)
        da_lep = real(delta_alpha(M_Z^2; alpha=α, fermions=SM_LEPTONS))
        da_had = real(delta_alpha(M_Z^2; alpha=α, fermions=SM_QUARKS))
        @test da_had > da_lep
    end

    @testset "Δα imaginary part (timelike)" begin
        # Above threshold 4m_e², Π̂ has a positive imaginary part (absorptive).
        # Since Δα = -Π̂ (Denner Eq. 3.10), Im(Δα) is NEGATIVE.
        da_mz = delta_alpha(M_Z^2; alpha=α)
        @test imag(da_mz) != 0.0  # timelike → complex
        @test imag(da_mz) < 0.0   # negative: Δα = -Π̂, absorptive Im(Π̂) > 0

        # Below all thresholds: purely real
        da_below = delta_alpha(1e-8; alpha=α)  # s << 4m_e²
        @test imag(da_below) ≈ 0.0 atol=1e-12
    end

    @testset "Δα energy dependence" begin
        # Δα increases with energy (QED is not asymptotically free)
        da_1 = real(delta_alpha(1.0; alpha=α))
        da_10 = real(delta_alpha(100.0; alpha=α))
        da_mz = real(delta_alpha(M_Z^2; alpha=α))
        @test da_1 < da_10 < da_mz

        # Spacelike q² < 0: Δα should be real (no thresholds crossed)
        da_sp = delta_alpha(-100.0; alpha=α)
        @test imag(da_sp) ≈ 0.0 atol=1e-10
        @test real(da_sp) > 0.0  # screening also for spacelike
    end

    # ----------------------------------------------------------------
    # Running α
    # ----------------------------------------------------------------
    @testset "Running α at M_Z" begin
        # PDG: α⁻¹(M_Z²) = 128.952 ± 0.014
        # Perturbative result will be close but not exact due to light quarks
        alpha_mz = running_alpha(M_Z^2; alpha=α)
        inv_alpha_mz = 1.0 / alpha_mz
        @test 127.0 < inv_alpha_mz < 131.0  # broad tolerance

        # Leptonic only: more precise
        alpha_mz_lep = running_alpha(M_Z^2; alpha=α, fermions=SM_LEPTONS)
        inv_lep = 1.0 / alpha_mz_lep
        @test 132.0 < inv_lep < 134.0  # α⁻¹ ≈ 137/(1-0.031) ≈ 133.4
    end

    @testset "Running α basic properties" begin
        # α(0) = α (low-energy limit)
        @test running_alpha(0.0; alpha=α) ≈ α rtol=1e-10

        # α increases with energy (screening)
        @test running_alpha(1.0; alpha=α) > α
        @test running_alpha(M_Z^2; alpha=α) > running_alpha(1.0; alpha=α)

        # Spacelike: running is real
        alpha_sp = running_alpha(-100.0; alpha=α)
        @test alpha_sp > α
        @test isfinite(alpha_sp)
    end

    @testset "Leading-log approximation" begin
        # For q² >> m² (single flavor, unit charge):
        # Δα ≈ (α/(3π)) × [ln(q²/m²) - 5/3]
        # Ref: Peskin & Schroeder, Ch. 7 (leading-log from schwinger.jl)
        m2 = 1.0
        q2 = 1e8  # q²/m² = 10⁸
        da_exact = real(delta_alpha(q2; alpha=α,
                        fermions=[(mass2=m2, charge2=1.0, nc=1)]))
        da_ll = α / (3π) * (log(q2 / m2) - 5 / 3)
        @test da_exact ≈ da_ll rtol=0.05  # 5% at s/m² = 10⁸
    end

    # ----------------------------------------------------------------
    # Improved Born cross-section
    # ----------------------------------------------------------------
    @testset "Improved Born σ(e⁺e⁻ → μ⁺μ⁻)" begin
        s = M_Z^2
        GEV2_TO_NB = 0.3894e6

        # Tree-level (in nb) and NLO without VP (convert to nb)
        sigma_tree_nb = 4π * α^2 / (3s) * GEV2_TO_NB
        sigma_nlo_nb = sigma_nlo_ee_mumu(s; alpha=α) * GEV2_TO_NB

        # Improved Born with VP (already in nb)
        sigma_imp = sigma_improved_ee_mumu(s; alpha=α)

        # VP increases effective α → improved > NLO
        @test sigma_imp > sigma_nlo_nb

        # Ratio: σ_improved/σ_tree ≈ (α(s)/α)² × (1 + Schwinger)
        alpha_s = running_alpha(s; alpha=α)
        ratio_expected = (alpha_s / α)^2 * (1 + schwinger_correction(; alpha=α))
        ratio_actual = sigma_imp / sigma_tree_nb
        @test ratio_actual ≈ ratio_expected rtol=1e-8

        # At low energy (s = 1 GeV²), VP effect is smaller
        sigma_imp_low = sigma_improved_ee_mumu(1.0; alpha=α)
        sigma_nlo_low_nb = sigma_nlo_ee_mumu(1.0; alpha=α) * GEV2_TO_NB
        @test sigma_imp_low / sigma_nlo_low_nb < sigma_imp / sigma_nlo_nb
    end

    @testset "σ improved at multiple energies" begin
        # Cross-section decreases as 1/s (up to logarithmic running)
        s_values = [10.0^2, 50.0^2, 91.2^2, 200.0^2]
        sigmas = [sigma_improved_ee_mumu(s; alpha=α) for s in s_values]

        # Monotonically decreasing with s
        for i in 1:length(sigmas)-1
            @test sigmas[i] > sigmas[i+1]
        end

        # All positive and finite
        @test all(isfinite, sigmas)
        @test all(s -> s > 0, sigmas)
    end
end
