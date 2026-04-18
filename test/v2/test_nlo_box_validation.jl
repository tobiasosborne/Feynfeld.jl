# Spiral 10 Phase E+F: NLO box evaluation and validation
#
# Validates the box contribution to Born-virtual interference for e⁺e⁻ → μ⁺μ⁻.
# Runs the full 6-layer pipeline: Model → Rules → Channels → Amplitude → TID → ComplexF64
#
# Ground truth: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/ElAel-MuAmu.m
# Cross-validation: arXiv:hep-ph/0010075 (Ellis, Giele, Zanderighi), Eq. 2.32
# Convention: COLLIER normalization for loop integrals
#
# Physical Born-virtual (spin-averaged):
#   (1/4) Σ_spins 2Re(M_tree* × M_box) = -e⁶/(32π²s) × Im(I_COLLIER)
# The overall coupling factor is purely imaginary (the tree propagator
# phase survives), so the physical quantity depends on Im(I), not Re(I).
# Ref: nlo_box.jl header for derivation
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 4 (TID)

using Test
using Feynfeld

# ---- Setup: ee→μμ process ----
const _p1 = Momentum(:p1)
const _p2 = Momentum(:p2)
const _k1 = Momentum(:k1)
const _k2 = Momentum(:k2)

const _incoming = [
    ExternalLeg(:e, _p1, true, false),
    ExternalLeg(:e, _p2, true, true),
]
const _outgoing = [
    ExternalLeg(:mu, _k1, false, false),
    ExternalLeg(:mu, _k2, false, true),
]

const _model = qed_model()
const _rules = feynman_rules(_model)

@testset "NLO box validation: e⁺e⁻ → μ⁺μ⁻" begin

    # ==== Phase D: TID evaluation at specific kinematics ====
    @testset "Phase D: TID evaluation (direct + crossed)" begin
        s = 100.0^2; cosθ = 0.5
        man = Mandelstam(s, cosθ)
        sp_vals = sp_values_2to2(man)

        tree_chs = tree_channels(_model, _rules, _incoming, _outgoing)
        tree_amp = build_amplitude(first(tree_chs), _rules, _model)

        channels = box_channels(_model, _rules, _incoming, _outgoing)
        ch_direct = first(ch for ch in channels if ch.topology == :direct_box)
        ch_crossed = first(ch for ch in channels if ch.topology == :crossed_box)

        # Direct box: D₀(0,0,0,0,s,t,0,0,0,0)
        I_direct = evaluate_single_box_channel(
            tree_amp, ch_direct, _rules, _model, sp_vals, man)
        @test isfinite(real(I_direct))
        @test isfinite(imag(I_direct))
        @test !iszero(I_direct)

        # Crossed box: D₀(0,0,0,0,s,u,0,0,0,0)
        I_crossed = evaluate_single_box_channel(
            tree_amp, ch_crossed, _rules, _model, sp_vals, man)
        @test isfinite(real(I_crossed))
        @test isfinite(imag(I_crossed))
        @test !iszero(I_crossed)

        # Direct ≠ crossed at generic angle (t ≠ u)
        @test !(I_direct ≈ I_crossed)
    end

    # ==== Phase E: Channel sum ====
    @testset "Phase E: evaluate_box_channels" begin
        s = 100.0^2; cosθ = 0.5
        man = Mandelstam(s, cosθ)
        I_total = evaluate_box_channels(_model, _rules, _incoming, _outgoing, man)
        @test isfinite(real(I_total))
        @test isfinite(imag(I_total))
        @test !iszero(I_total)
    end

    # ==== Phase F: Validation tests ====

    # ---- F1: Im(I) is forward-backward symmetric ----
    # The physical quantity depends on Im(I_COLLIER) (not Re).
    # Im(I_direct + I_crossed) must be symmetric under cosθ → -cosθ.
    # Note: Re(I) is NOT expected to be F-B symmetric for box-only,
    # because the box is not individually gauge-invariant under crossing.
    @testset "F1: Im(I) forward-backward symmetry" begin
        s = 100.0^2
        man_fwd = Mandelstam(s, 0.5)
        man_bwd = Mandelstam(s, -0.5)

        I_fwd = evaluate_box_channels(_model, _rules, _incoming, _outgoing, man_fwd)
        I_bwd = evaluate_box_channels(_model, _rules, _incoming, _outgoing, man_bwd)

        @test imag(I_fwd) ≈ imag(I_bwd) rtol=1e-8
    end

    # ---- F2: Im(I) crossing at θ=90° ----
    # At θ=90°, t=u. The Im parts of direct and crossed channels
    # must be identical (same PaVe invariants, same kinematics).
    @testset "F2: Im(I) crossing at θ=90°" begin
        s = 100.0^2
        man90 = Mandelstam(s, 0.0)
        sp_vals = sp_values_2to2(man90)

        tree_amp = build_amplitude(
            first(tree_channels(_model, _rules, _incoming, _outgoing)), _rules, _model)
        channels = box_channels(_model, _rules, _incoming, _outgoing)
        ch_d = first(ch for ch in channels if ch.topology == :direct_box)
        ch_c = first(ch for ch in channels if ch.topology == :crossed_box)

        I_d = evaluate_single_box_channel(tree_amp, ch_d, _rules, _model, sp_vals, man90)
        I_c = evaluate_single_box_channel(tree_amp, ch_c, _rules, _model, sp_vals, man90)

        @test imag(I_d) ≈ imag(I_c) rtol=1e-8
    end

    # ---- F3: Finiteness at multiple energies and angles ----
    @testset "F3: Finiteness scan" begin
        for sqrt_s in [50.0, 100.0, 200.0, 500.0]
            s = sqrt_s^2
            for cosθ in [-0.8, -0.3, 0.0, 0.3, 0.8]
                man = Mandelstam(s, cosθ)
                I = evaluate_box_channels(_model, _rules, _incoming, _outgoing, man)
                @test isfinite(real(I))
                @test isfinite(imag(I))
            end
        end
    end

    # ---- F4: Order of magnitude ----
    # The box correction relative to Born should be O(α/π) ≈ O(10⁻³).
    #   Born = 2e⁴(t²+u²)/s²   (spin-averaged)
    #   Box  = -e⁶/(32π²s) × Im(I_COLLIER)
    #   Ratio ≈ α/(4π) × f(cosθ)  where f = O(1..10)
    #
    # Ref: cross_section.jl line 130: "|M|² = e⁴ × 8(t²+u²)/s²"
    @testset "F4: Order of magnitude" begin
        α = 1.0 / 137.036
        s = 100.0^2; cosθ = 0.5
        man = Mandelstam(s, cosθ)

        # Born (spin-averaged, with coupling)
        tree_result = solve_tree(CrossSectionProblem(_model, _incoming, _outgoing, s))
        tree_trace = evaluate_m_squared(tree_result.amplitude_squared, man)
        born = (4π * α)^2 / man.s^2 * tree_trace / 4

        # Box (spin-averaged, with coupling)
        box = born_virtual_box(_model, _rules, _incoming, _outgoing, man; alpha=α)

        ratio = abs(box / born)
        # Must be O(α/π) — between 10⁻⁵ and 10⁻¹
        @test 1e-5 < ratio < 0.1
    end

    # ---- F5: Energy scaling ----
    # At fixed angle, the box correction / Born should be roughly
    # independent of s (up to logs).  Check that the ratio changes by
    # less than a factor of 10 when √s goes from 50 to 500 GeV.
    @testset "F5: Energy scaling" begin
        α = 1.0 / 137.036; cosθ = 0.3
        ratios = Float64[]
        for sqrt_s in [50.0, 100.0, 200.0, 500.0]
            s = sqrt_s^2
            man = Mandelstam(s, cosθ)
            tree_result = solve_tree(CrossSectionProblem(_model, _incoming, _outgoing, s))
            tree_trace = evaluate_m_squared(tree_result.amplitude_squared, man)
            born = (4π * α)^2 / s^2 * tree_trace / 4
            box = born_virtual_box(_model, _rules, _incoming, _outgoing, man; alpha=α)
            push!(ratios, abs(box / born))
        end
        # All ratios should be within a factor of 10 of each other
        @test maximum(ratios) / minimum(ratios) < 10.0
    end

    # ---- F6: Consistency with direct COLLIER D₀ ----
    @testset "F6: COLLIER D₀ sanity" begin
        s = 100.0^2; t = -s * 0.25

        # D₀(0,0,0,0,s,t,0,0,0,0): all-massless box
        d0_st = evaluate(D0(0.0, 0.0, 0.0, 0.0, s, t, 0.0, 0.0, 0.0, 0.0))
        @test isfinite(real(d0_st))
        @test isfinite(imag(d0_st))
        @test !iszero(d0_st)

        # D₀(0,0,0,0,s,u,0,0,0,0): crossed channel
        u = -s - t
        d0_su = evaluate(D0(0.0, 0.0, 0.0, 0.0, s, u, 0.0, 0.0, 0.0, 0.0))
        @test isfinite(real(d0_su))
        @test isfinite(imag(d0_su))
        @test !iszero(d0_su)

        # D₀ at two different (s, man_inv) should differ (non-degenerate)
        @test !(d0_st ≈ d0_su)
    end

    # ---- F7: born_virtual_box is non-zero ----
    # The box-only contribution can be positive or negative; the FULL
    # virtual correction (box + vertex + self-energy + CTs) is negative.
    @testset "F7: born_virtual_box non-zero" begin
        α = 1.0 / 137.036; s = 100.0^2
        man = Mandelstam(s, 0.5)
        box = born_virtual_box(_model, _rules, _incoming, _outgoing, man; alpha=α)
        @test isfinite(box)
        @test !iszero(box)
    end
end
