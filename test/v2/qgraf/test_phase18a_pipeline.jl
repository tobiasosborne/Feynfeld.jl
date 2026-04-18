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

    @testset "ee→μμ tree (massless): pipeline == solve_tree" begin
        in_e   = ExternalLeg(:e, Momentum(:p1), true,  false)
        in_eb  = ExternalLeg(:e, Momentum(:p2), true,  true)
        out_mu = ExternalLeg(:mu, Momentum(:k1), false, false)
        out_mb = ExternalLeg(:mu, Momentum(:k2), false, true)
        prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                                    [in_e, in_eb], [out_mu, out_mb], 10.0)

        handbuilt = solve_tree(prob)
        pipeline  = solve_tree_pipeline(prob)

        @test pipeline.amplitude_squared == handbuilt.amplitude_squared
    end

end
