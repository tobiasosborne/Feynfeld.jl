#  Phase 14a: inline-simple xg-pattern filters.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                              has_no_selfloop, has_no_diloop, has_no_parallel

@testset "Phase 14a: inline xg filters" begin

    "Build a TopoState with caller-controlled xg entries."
    function _state(n_ext, counts, mrho, nloop)
        TopoState(Partition(Int8(n_ext), Int8.(counts), Int8(mrho), Int8(nloop)))
    end

    @testset "has_no_selfloop" begin
        s = _state(0, [2], 3, 2)
        # Empty xg ⇒ no self-loops.
        @test has_no_selfloop(s) == true
        s.xg[1, 1] = Int8(2)   # one self-loop on vertex 1
        @test has_no_selfloop(s) == false
        # Externals' diagonal doesn't count (loop only over internals).
        s2 = _state(2, [1], 3, 0)
        s2.xg[1, 1] = Int8(2)  # external — irrelevant
        @test has_no_selfloop(s2) == true
    end

    @testset "has_no_diloop" begin
        s = _state(0, [2], 3, 2)
        @test has_no_diloop(s) == true
        s.xg[1, 2] = Int8(2)   # double edge between 1 and 2
        @test has_no_diloop(s) == false
        s.xg[1, 2] = Int8(1)
        @test has_no_diloop(s) == true
    end

    @testset "has_no_parallel" begin
        # Same as nodiloop plus xg[i,i] ≤ 3 for internals.
        s = _state(0, [2], 3, 2)
        @test has_no_parallel(s) == true
        s.xg[1, 1] = Int8(4)   # quadruple self-loop pair (forbidden)
        @test has_no_parallel(s) == false
        s.xg[1, 1] = Int8(2)
        @test has_no_parallel(s) == true
        # Diloop also rejects.
        s.xg[1, 2] = Int8(2)
        @test has_no_parallel(s) == false
    end

end
