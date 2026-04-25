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
# Session 32 fix (bead feynfeld-vjw9 closed): the orbit-rep dedup bug was
# in the ps1 action used by `is_emission_canonical` / `emission_stabilizer`
# / `same_emission_orbit`.  All three used the LEFT action
# `(g·ps1)[i] = g[ps1[i]]` (relabels physical legs), but `_diagram_sig`
# uses the right action on vertex relabeling.  Switching ps1 to the right
# action `(g·ps1)[i] = ps1[g⁻¹[i]]` (relabels SLOTS — the physically
# correct convention since topology autos describe combinatorial symmetry,
# not physical leg reidentification) made all three agree, and
# `is_emission_canonical` now accepts exactly `count_diagrams_qg21(...)`
# emissions per process — verified for Bhabha (2), ee→μμ (1).
#
# With 2 canonical reps for Bhabha, the interference loop in
# `spin_sum_interference` closes cleanly across the s/t reps because their
# bar_mom assignments tile a 4-cycle (p1 → k1 → k2 → p2 → p1), so no
# bar_mom canonicalisation (Option B / bead `feynfeld-rj1l`) is needed
# at this scope.  Both former @test_broken assertions now pass.

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
        @test result.n_emissions == 2
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

        @test result.amplitude_squared == expected_expanded
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
