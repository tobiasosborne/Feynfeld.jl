#  Phase 18a-8: pipeline integration solve_tree_pipeline.
#
#  Smoke tests that solve_tree_pipeline drives the qg21 emission
#  stream through emission_to_amplitude into the existing spin-sum /
#  contract / expand machinery and returns a non-trivial AlgSum.
#  Symbolic equality with the existing solve_tree path is the
#  acceptance test for Phase 18a-9.

using Test

using Feynfeld

@testset "Phase 18a-8: solve_tree_pipeline" begin

    @testset "ee→μμ tree (massless): pipeline returns a non-zero AlgSum" begin
        in_e   = ExternalLeg(:e, Momentum(:p1), true,  false)
        in_eb  = ExternalLeg(:e, Momentum(:p2), true,  true)
        out_mu = ExternalLeg(:mu, Momentum(:k1), false, false)
        out_mb = ExternalLeg(:mu, Momentum(:k2), false, true)
        prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                                    [in_e, in_eb], [out_mu, out_mb], 10.0)
        result = solve_tree_pipeline(prob)
        @test result.amplitude_squared isa AlgSum
        @test !iszero(result.amplitude_squared)
        @test result.n_emissions >= 1
    end

end
