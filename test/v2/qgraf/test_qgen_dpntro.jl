#  Phase 12a: dpntro lookup table for qgen field assignment.
#
#  qgraf's dpntro is an arena-indexed lookup: dpntro(degree) → first-particle
#  index → list of vertex rules (sorted field tuples).  The Julia equivalent
#  is a nested Dict.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13889
#    "vfo(vv) = stib(stib(dpntro(0)+vdeg(vv))+pmap(vv,1))"
#  This is a two-level lookup: first by degree, then by leading field index.
#
#  For our port we use Dict{Int, Vector{Vector{Symbol}}} (degree → sorted
#  list of expanded vertex rules); a future optimisation can split into
#  the nested form when we measure a hot path.

using Test
using Feynfeld
using Feynfeld: _expand_model_for_diagen
using Feynfeld.QgrafPort: build_dpntro

@testset "Phase 12a: dpntro lookup table" begin

    @testset "phi3 model: single 3-arity vertex" begin
        m = phi3_model()
        exp = _expand_model_for_diagen(m)
        dp = build_dpntro(exp.vertex_rules)
        @test haskey(dp, 3)
        @test length(dp[3]) == 1
        @test dp[3][1] == sort([:phi, :phi, :phi])
    end

    @testset "QED model: 2-gen has 2 cubic vertices" begin
        # qed_model() = 2-gen (e + μ).  Two qqg vertices: eeγ and μμγ.
        m = qed_model()
        exp = _expand_model_for_diagen(m)
        dp = build_dpntro(exp.vertex_rules)
        @test haskey(dp, 3)
        # _expand_vertex maps (:e,:e,:gamma) → sorted(:e, :e_bar, :gamma)
        # and (:mu,:mu,:gamma) → sorted(:mu, :mu_bar, :gamma).
        @test length(dp[3]) == 2
    end

    @testset "QCD model: cubic + quartic + ghost" begin
        m = qcd_model()
        exp = _expand_model_for_diagen(m)
        dp = build_dpntro(exp.vertex_rules)
        # qqg cubic, ggg cubic, ghost-ghost-g cubic → 3 cubic rules
        # gggg quartic → 1 quartic rule
        @test haskey(dp, 3)
        @test haskey(dp, 4)
        @test length(dp[3]) == 3
        @test length(dp[4]) == 1
        @test dp[4][1] == [:g, :g, :g, :g]
    end

    @testset "All rules in dp[d] have length d" begin
        for m in [phi3_model(), qed_model(), qcd_model()]
            exp = _expand_model_for_diagen(m)
            dp = build_dpntro(exp.vertex_rules)
            for (d, rules) in dp
                for rule in rules
                    @test length(rule) == d
                end
            end
        end
    end

    @testset "All rules in dp[d] are sorted" begin
        for m in [phi3_model(), qed_model(), qcd_model()]
            exp = _expand_model_for_diagen(m)
            dp = build_dpntro(exp.vertex_rules)
            for (_, rules) in dp
                for rule in rules
                    @test issorted(rule)
                end
            end
        end
    end

end
