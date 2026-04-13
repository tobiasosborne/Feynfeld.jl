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

    # ── Phase 7b+7c: dsum self-loop loop + row-by-row off-diag fill ─────

    @testset "Phase 7c: phi3 φφ→φφ tree — exactly 1 topology" begin
        # 4 ext, 2 deg-3 internals, nloop=0.
        # Hand trace through qg21:
        #   Step B emits 5 (xc, xn); 3 rejected by rho<n filter.
        #   Surviving: xn=[3,1] and xn=[2,2].
        #   Step C: xn=[3,1] fails row-fill (vertex-6 needs 2 edges but bond=1).
        #          xn=[2,2] yields the unique tree topology xg[5,6]=1.
        # qgraf "3 diagrams" for φφ→φφ tree comes from external-leg labelling
        # (qg10) producing 3 channels (s, t, u) atop this single topology.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = _count_emissions(s)
        @test length(emits) == 1
        # The single topology: vertex 5 ↔ vertex 6 with 1 internal edge.
        @test emits[1][5, 6] == Int8(1)
    end

    @testset "Phase 7b: dsum>0 in vacuum 2-loop dumbbell" begin
        # 0 ext, 2 deg-3 internals, nloop=2.  All edges between vertices 1,2.
        # vdeg=[3,3], P=3, V=2, L=2 ✓
        # Topologies:
        #   - 3 parallel edges between v1-v2 (xg[1,2]=3, xg[1,1]=xg[2,2]=0)
        #   - 1 edge v1-v2 + each has 1 self-loop (xg[1,2]=1, xg[1,1]=2, xg[2,2]=2)
        # qg21 raw count: 2 topologies.
        p = Partition(Int8(0), Int8[2], Int8(3), Int8(2))
        s = TopoState(p)
        emits = _count_emissions(s)
        @test length(emits) == 2
    end

    @testset "Phase 7c: QED ee→μμ-shape tree — 1 topology" begin
        # Same partition shape as phi3 φφ→φφ tree (4 ext, 2 deg-3, nloop=0).
        # Confirms the shape-equivalence: at qg21 level the topology count
        # is the same; the FIELD ASSIGNMENT in qgen distinguishes ee→μμ
        # from φφ→φφ.  Field-level diagram count for ee→μμ tree = 1
        # (matches existing test_ee_mumu_x.jl).
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = _count_emissions(s)
        @test length(emits) == 1
    end

end
