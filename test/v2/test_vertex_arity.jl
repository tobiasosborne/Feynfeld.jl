#  Test: variable-arity VertexRule.
#
#  qgraf models declare vertices of arbitrary arity:
#     [ phi, phi, phi ]                   3-vertex
#     [ gluon, gluon, gluon, gluon ]      4-vertex
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd verbatim
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/phi3 verbatim
#
#  Phase 2 of the Strategy C qgraf port: relax VertexRule.fields from
#  NTuple{3,Symbol} to Tuple{Vararg{Symbol}} so QCD's 4-gluon vertex,
#  ghost-gluon, and any future higher-arity coupling can be expressed.

using Test
using Feynfeld
using Feynfeld: VertexRule, FeynmanRules, feynman_rules

@testset "VertexRule variable arity" begin

    @testset "3-arity (regression)" begin
        v = VertexRule((:phi, :phi, :phi), :g_phi3)
        @test v.fields == (:phi, :phi, :phi)
        @test length(v.fields) == 3
        @test v.coupling === :g_phi3
    end

    @testset "4-arity (qgraf 4-gluon)" begin
        # Ref: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd
        # "[ gluon, gluon, gluon, gluon ]"
        v = VertexRule((:g, :g, :g, :g), :g_s)
        @test v.fields == (:g, :g, :g, :g)
        @test length(v.fields) == 4
        @test v.coupling === :g_s
    end

    @testset "5-arity (forward compat)" begin
        v = VertexRule((:a, :b, :c, :d, :e), :coup)
        @test length(v.fields) == 5
    end

end

@testset "FeynmanRules.vertices variable-arity Dict" begin
    # Build a QCD-shaped FeynmanRules by hand: 3-arity qqg + 4-arity gggg.
    @testset "mixed-arity Dict construction" begin
        verts = Dict{Tuple{Vararg{Symbol}}, VertexRule}()
        verts[(:q, :q, :g)] = VertexRule((:q, :q, :g), :g_s)
        verts[(:g, :g, :g)] = VertexRule((:g, :g, :g), :g_s)
        verts[(:g, :g, :g, :g)] = VertexRule((:g, :g, :g, :g), :g_s)
        @test length(verts) == 3
        @test haskey(verts, (:g, :g, :g, :g))
        @test verts[(:g, :g, :g, :g)].fields == (:g, :g, :g, :g)
    end
end
