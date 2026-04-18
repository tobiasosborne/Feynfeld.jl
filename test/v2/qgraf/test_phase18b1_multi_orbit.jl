# Phase 18b-1 (bead feynfeld-ewgw) + Phase 18b-1a (bead feynfeld-vjw9):
# multi-orbit Burnside summation validated on Bhabha.
#
# Ground truth:
#   - qgraf count_diagrams_qg21(QED ee→ee tree) = 2 (s-channel annihilation
#     + t-channel exchange).  Session 22 Phase 17a VERDICT table (Burnside
#     strategy, 9/10).
#   - Handbuilt trace-only |M|² = T_tt + T_ss − 2·T_int with fermion signs
#     (+, −) for (t, s) from the identical-fermion Fermi-Dirac relative
#     minus.  Reference implementation: test/v2/test_bhabha.jl,
#     FeynCalc Examples/QED/Tree/Mathematica/ElAel-ElAel.m:98-99.
#
# solve_tree_pipeline must produce n_emissions == 2 (one canonical rep per
# orbit) and amplitude_squared symbolically equal to the handbuilt sum.
#
# Session 29 status: the `@test_broken` assertions below are the RED tests
# for bead feynfeld-vjw9.  The Session 27 `is_emission_canonical` filter
# under-counts (Strategy-C bug, HANDOFF Session 22): canonical lex-smallest
# rep of one Bhabha orbit is qgen-invalid, so that orbit yields 0 instead
# of 1 emission.  Session 29 tried Burnside weighting (every emission ×
# |Stab|/|G|) — arithmetic works, but `_find_line_by_bar_mom` inside
# spin_sum_interference errors on cross-bundle automorphic relabelings
# ("No line with bar momentum ..."). A new helper `same_emission_orbit`
# (audition.jl) partitions orbits correctly for ee→μμ (1 orbit) but splits
# Bhabha into 4 orbits instead of 2 — further auto-action calibration is
# required. Next attempt: either fix interference label matching (bead's
# Option B) to make Burnside-all work, or finish same_emission_orbit (bead's
# Option C). Either unlocks these @test_broken assertions.

using Test
using Feynfeld
using Feynfeld.QgrafPort: count_diagrams_qg21

@testset "Phase 18b-1: multi-orbit Bhabha" begin
    p1 = Momentum(:p1); p2 = Momentum(:p2); k1 = Momentum(:k1); k2 = Momentum(:k2)

    @testset "qgraf reference: Bhabha has 2 orbits" begin
        @test count_diagrams_qg21(qed_model(m_e=:zero),
                                    [:e, :e_bar], [:e, :e_bar]) == 2
    end

    in_e   = ExternalLeg(:e,     p1, true,  false)
    in_eb  = ExternalLeg(:e,     p2, true,  true)
    out_e  = ExternalLeg(:e,     k1, false, false)
    out_eb = ExternalLeg(:e,     k2, false, true)
    prob = CrossSectionProblem(qed_model(m_e=:zero),
                                [in_e, in_eb], [out_e, out_eb], 10.0)

    result = solve_tree_pipeline(prob)

    @testset "solve_tree_pipeline produces one rep per orbit" begin
        @test_broken result.n_emissions == 2
    end

    @testset "trace-only |M|² equals handbuilt T_tt + T_ss − 2 T_int" begin
        # Build the handbuilt trace-only reference using the same Dirac-trace
        # machinery the pipeline uses downstream.  Ref: test/v2/test_bhabha.jl.
        Tr_t1 = dirac_trace(DiracGamma[GS(p1), GAD(:alpha_), GS(k1), GAD(:alpha)])
        Tr_t2 = dirac_trace(DiracGamma[GS(k2), GAD(:alpha_), GS(p2), GAD(:alpha)])
        T_tt  = Tr_t1 * Tr_t2

        Tr_s1 = dirac_trace(DiracGamma[GS(p1), GAD(:beta_),  GS(p2), GAD(:beta)])
        Tr_s2 = dirac_trace(DiracGamma[GS(k2), GAD(:beta_),  GS(k1), GAD(:beta)])
        T_ss  = Tr_s1 * Tr_s2

        # Cross-line 8-gamma closed trace (M_s* × M_t).
        T_int = dirac_trace(DiracGamma[
            GS(p1), GAD(:beta_),  GS(p2), GAD(:alpha),
            GS(k2), GAD(:beta_),  GS(k1), GAD(:alpha),
        ])

        expected_raw      = T_tt + T_ss - 2 * T_int
        expected_expanded = expand_scalar_product(contract(expected_raw))

        @test_broken result.amplitude_squared == expected_expanded
    end

    @testset "Phase 18a regression: ee→μμ still single-orbit" begin
        in_e2   = ExternalLeg(:e,  Momentum(:p1), true,  false)
        in_eb2  = ExternalLeg(:e,  Momentum(:p2), true,  true)
        out_mu  = ExternalLeg(:mu, Momentum(:k1), false, false)
        out_mub = ExternalLeg(:mu, Momentum(:k2), false, true)
        prob_mumu = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                                          [in_e2, in_eb2], [out_mu, out_mub], 10.0)
        r = solve_tree_pipeline(prob_mumu)
        @test r.n_emissions == 1
    end
end
