# Single-process test runner for Feynfeld.
# Loads Feynfeld once, then runs all test files in the same Julia process.
#
# Usage:
#   julia --project=. test/v2/runtests.jl
#   julia --project=. -e 'using Pkg; Pkg.test()'  (via test/runtests.jl forwarder)

using Test

using Feynfeld

@testset "Feynfeld" begin
    # ---- Core pipeline / algebra / integrals ----
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
    # ---- Move 1.3: formerly missing from orchestration ----
    include("test_diagram_gen.jl")
    include("test_vertex_arity.jl")
    include("test_qcd_4gluon.jl")
    include("test_qcd_ghost.jl")
    include("test_ee_ww_grozin.jl")
    # ---- MUnit translations (FeynCalc port) ----
    for f in sort(readdir(joinpath(@__DIR__, "munit")))
        startswith(f, "test_") && endswith(f, ".jl") || continue
        include(joinpath("munit", f))
    end
    # ---- qgraf Strategy C port integration tests ----
    for f in sort(readdir(joinpath(@__DIR__, "qgraf")))
        startswith(f, "test_") && endswith(f, ".jl") || continue
        include(joinpath("qgraf", f))
    end
end
