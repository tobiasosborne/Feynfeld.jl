#  Phase 14c: qumvi family — vertex-cut filters (nosnail, onevi).
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:3777-3881.
#  qgsig (nosigma) and qcyc (cycli) depend on flow[][] — deferred to Phase 16.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                              has_no_snail, is_one_vi

@testset "Phase 14c: qumvi vertex filters" begin

    @testset "has_no_snail: tree-level always passes" begin
        # qg21:3789-3791 — nloop==0 → trivially OK (no snails possible).
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        @test has_no_snail(s) == true
    end

    @testset "has_no_snail: 1-ext 1-loop self-energy: snail IS the diagram" begin
        # qg21:3792-3797 — exception: rho(-1)==1 AND nloop==1 → don't reject.
        # Construct a 1-loop self-energy: 1 ext + 1 internal with self-loop.
        # Need V valid: 1 + 1 = 2, half-edges = 1 + 3 = 4 → P=2, V=2, L=1 ✓.
        p = Partition(Int8(1), Int8[1], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[2, 2] = Int8(2)        # self-loop on internal vertex 2
        @test has_no_snail(s) == true
    end

    @testset "has_no_snail: 2-ext 1-loop with snail → reject" begin
        # 2 ext + 2 internals + 1L: half-edges = 2 + 6 = 8, P=4, V=4, L=1.
        # If vertex 4 has self-loop and is connected only via single bridge:
        # rho(-1)==2 (NOT 1), nloop==1 → exception doesn't apply.
        p = Partition(Int8(2), Int8[2], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[3, 3] = Int8(2)        # self-loop on internal vertex 3
        @test has_no_snail(s) == false
    end

    @testset "is_one_vi: 1 internal vertex → vacuously 1VI" begin
        # qg21:3806-3808 — n - 1 - rhop1 <= 0.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        @test is_one_vi(s) == true
    end

    @testset "is_one_vi: phi3 φφ→φφ tree has 2 internals — neither is a cut" begin
        # 4 ext + 2 internals (5,6); single internal edge xg[5,6].
        # Removing vertex 5: vertex 6 is alone among internals; OK (no cut).
        # Removing vertex 6: vertex 5 alone; OK.
        # Both internals isolated when other removed — no cut.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        s.xn[5] = Int8(2); s.xn[6] = Int8(2)
        s.xg[1, 5] = Int8(1); s.xg[2, 5] = Int8(1)
        s.xg[3, 6] = Int8(1); s.xg[4, 6] = Int8(1)
        s.xg[5, 6] = Int8(1)
        @test is_one_vi(s) == true
    end

    @testset "is_one_vi: middle vertex on a 3-internal chain IS a cut" begin
        # 4 ext + 3 deg-3 internals + 1L? half-edges = 4+9=13 odd ✗
        # 2 ext + 3 deg-3 internals + 1L: 2+9=11 odd ✗
        # 4 ext + 2 deg-3 + 1 deg-4 internals + 1L: 4+6+4=14, P=7, V=7, L=1 ✓.
        # Layout: ext{1,2}→v5; ext{3,4}→v6; v5—v7—v6 chain via the deg-4 vertex
        # in the middle.  Plus an internal-edge cycle: v5—v7, v6—v7 (2 edges
        # to v7), and v7 has a self-loop (xg[7,7]=2).  Needs to satisfy degrees.
        # Simpler: chain v5—v7—v6, v7 has degree 4; needs (5—7), (6—7), and 2
        # self-loop half-edges (xg[7,7]=2).
        # vdeg: v5=3 (xn=2 + xg[5,7]=1), v6=3 (xn=2 + xg[6,7]=1), v7=4 (xg[5,7]+xg[6,7]+xg[7,7]=1+1+2=4) ✓.
        p = Partition(Int8(4), Int8[2, 1], Int8(3), Int8(1))
        s = TopoState(p)
        s.xn[5] = Int8(2); s.xn[6] = Int8(2); s.xn[7] = Int8(0)
        s.xg[1, 5] = Int8(1); s.xg[2, 5] = Int8(1)
        s.xg[3, 6] = Int8(1); s.xg[4, 6] = Int8(1)
        s.xg[5, 7] = Int8(1); s.xg[6, 7] = Int8(1)
        s.xg[7, 7] = Int8(2)
        # Removing v7: v5 cannot reach v6 → 1VI fails.
        @test is_one_vi(s) == false
    end

end
