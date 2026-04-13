#  Phase 4: ghost field + ghost-gluon vertex in qcd_model.
#
#  Faddeev-Popov ghosts are anticommuting (fermion statistics) but Lorentz
#  scalars; they cancel longitudinal-gluon contributions in covariant gauges.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd
#     [ ghost, antighost, -1 ]            ! propagator (fermion statistics)
#     [ antighost, ghost, gluon ]         ! ghost-gluon vertex
#
#  For diagram-counting purposes we model the ghost as Field{Fermion}; the
#  Lorentz-scalar treatment is a Layer-4 concern (TBD).

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX: model_fields, feynman_rules

@testset "QCD ghost field + ghost-gluon vertex" begin

    qcd = qcd_model()

    @testset "ghost field present in model" begin
        names = Symbol[f.name for f in model_fields(qcd)]
        @test :ghost in names
    end

    @testset "ghost-ghost-gluon vertex registered" begin
        rules = feynman_rules(qcd)
        @test haskey(rules.vertices, (:ghost, :ghost, :g))
    end

    @testset "ghost → ghost 1L onepi (golden master)" begin
        # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/qcd/
        #         ghost_TO_ghost_1L_onepi.out → 1 diagram
        @test count_diagrams(qcd, [:ghost], [:ghost]; loops=1, onepi=true) == 1
    end

end
