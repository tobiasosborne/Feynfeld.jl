#  Phase 7: qg21 Step C — main topology generation.
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12659-13150.
#
#  Composed of three sub-phases:
#    7a — external edge placement + dsum=0 + trivial single-internal emission
#    7b — dsum self-loop enumeration + xg-diagonal backtracking
#    7c — row-by-row off-diagonal xg fill + canonical cross-row check
#
#  Tests are added incrementally as each sub-phase lands.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, qg21_enumerate!

@testset "qg21 Step C" begin

    "Test helper: count emissions and snapshot adjacency matrices."
    function _count_emissions(state)
        emits = Matrix{Int8}[]
        qg21_enumerate!(state) do s
            n = Int(s.n)
            push!(emits, copy(s.xg[1:n, 1:n]))
        end
        emits
    end

    # ── Phase 7a: trivial single-internal cases ──────────────────────────

    @testset "Phase 7a: phi3 tree φ→φφ — exactly 1 topology" begin
        # 3 externals, 1 deg-3 internal, nloop=0.
        # Step B emits xc=[0,3] xn=[0,0,0,3] (xc=[2,1] rejected by rho<n filter).
        # Step C: place xg(j,4)=1 for j=1..3.  No row fill needed (n_int=1).
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        emits = _count_emissions(s)
        @test length(emits) == 1
        # The unique adjacency: vertex 4 connects to all three externals.
        adj = emits[1]
        @test adj[1, 4] == Int8(1)
        @test adj[2, 4] == Int8(1)
        @test adj[3, 4] == Int8(1)
        @test adj[4, 4] == Int8(0)
    end

    @testset "Phase 7a: rejects xc(1)>0 when rho<n (qg21:12660-12664)" begin
        # Use the same partition; if Step B emits xc=[2,1] xn=[0,0,0,1] the
        # filter must reject (rho(-1)=3 < n=4 AND xc(1)=2 > 0).  The total
        # emission count from Step C must therefore equal 1 even though
        # Step B alone would yield 2 emissions.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        @test length(_count_emissions(s)) == 1
    end

end
