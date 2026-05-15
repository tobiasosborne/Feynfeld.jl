#  Phase 18a-9: symbolic-equality validation.
#
#  For each test process, build |M|² two ways:
#    (a) the existing channels.jl / amplitude.jl / cross_section.solve_tree
#        hand-built path
#    (b) the new qg21 pipeline via solve_tree_pipeline
#  and assert AlgSum equality after spin-sum / contract / expand_sp.
#
#  Tobias Rule 7: symbolic only — no numerical spot-checks. Exact
#  Rational coefficients throughout.

using Test

using Feynfeld

@testset "Phase 18a-9: pipeline ≡ handbuilt symbolic equality" begin

    @testset "ee→μμ tree (massless): pipeline == solve_tree (kinematic)" begin
        # Bridge test (Phase 18b-7, feynfeld-5d1k): with coupling now flowing
        # through the pipeline (e^4 for ee→μμ), the bridge invariant becomes
        # `pipeline_kinematic(pipeline) == solve_tree` — `solve_tree` doesn't
        # carry the coupling factor, so we strip it from the pipeline result
        # for the kinematic comparison and assert the coupling separately.
        # solve_tree is marked for retirement once Phase 18b is complete.
        in_e   = ExternalLeg(:e, Momentum(:p1), true,  false)
        in_eb  = ExternalLeg(:e, Momentum(:p2), true,  true)
        out_mu = ExternalLeg(:mu, Momentum(:k1), false, false)
        out_mb = ExternalLeg(:mu, Momentum(:k2), false, true)
        prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                                    [in_e, in_eb], [out_mu, out_mb], 10.0)

        handbuilt = solve_tree(prob)
        pipeline  = solve_tree_pipeline(prob)

        @test pipeline_kinematic(pipeline).amplitude_squared ==
              handbuilt.amplitude_squared
        # Coupling assertion: every term carries e^4 (two vertices per fermion
        # line, two fermion lines in |M|² = M × M*). Ref: P&S §4.8.
        e4 = coupling_alg(:e, 4)
        @test pipeline.amplitude_squared ==
              e4 * pipeline_kinematic(pipeline).amplitude_squared
    end

end
