# Single-process test runner for Feynfeld v2.
# Loads FeynfeldX once, then runs all test files in the same Julia process.
# Individual test files still work standalone via @isdefined guard.
#
# Usage: julia --project=. test/v2/runtests.jl

using Test

include(joinpath(@__DIR__, "..", "..", "src", "v2", "FeynfeldX.jl"))
using .FeynfeldX

@testset "Feynfeld v2" begin
    include("test_coeff.jl")
    include("test_colour.jl")
    include("test_ee_mumu_x.jl")
    include("test_self_energy.jl")
    include("test_vertical.jl")
    include("test_pave.jl")
    include("test_schwinger.jl")
    include("test_compton.jl")
    include("test_munit_batch1.jl")
    include("test_munit_batch2.jl")
    include("test_bhabha.jl")
    include("test_qqbar_gg.jl")
    include("test_self_energy_1loop.jl")
    include("test_d0.jl")
    include("test_running_alpha.jl")
    include("test_ee_ww.jl")
    include("test_pipeline.jl")
    include("test_vertex_g2.jl")
    include("test_box_ee_mumu.jl")
    include("test_nlo_box_validation.jl")
end
