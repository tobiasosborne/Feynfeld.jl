#  Phase 18a-2: half-edge labelling (compute_amap).
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12133-12158
#  External back-writes: f08:12342-12344
#  Pairing invariant: amap[v, s] == amap[vmap[v, s], lmap[v, s]]
#
#  amap[v, s] = global half-edge id (1..nleg for externals, nleg+1..nli for
#  internals).  Both half-edges of any single edge carry the same id.
#  Edge-id ordering follows ege[i,j] = lex (i,j) for internal pairs.

using Test

include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels, compute_amap

@testset "Phase 18a-2: half-edge labelling (amap)" begin

    @testset "φ→φφ tree: externals + back-writes" begin
        # 3 ext (1,2,3) → 1 int (4).  No internal half-edges.
        # Expected: externals get amap[i,1] = i; their internal-vertex peer
        # slots get the same id back via the f08:12342-12344 back-write.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 4] = 1
        s.xg[2, 4] = 1
        s.xg[3, 4] = 1
        s.xn[4]    = 3
        labels = compute_qg10_labels(s)

        amap = compute_amap(s, labels)

        @test size(amap) == (4, MAX_V)
        @test amap[1, 1] == 1
        @test amap[2, 1] == 2
        @test amap[3, 1] == 3
        # Back-writes to internal vertex 4
        @test amap[4, 1] == 1
        @test amap[4, 2] == 2
        @test amap[4, 3] == 3
    end

    @testset "ee→μμ s-channel: internal edge has id nleg+1" begin
        # 4 ext, 2 int (5,6). One internal edge between 5 and 6.
        # ege[5,6] = rhop1 = nleg+1 = 5.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 5] = 1
        s.xg[2, 5] = 1
        s.xg[3, 6] = 1
        s.xg[4, 6] = 1
        s.xg[5, 6] = 1
        s.xn[5]    = 2
        s.xn[6]    = 2
        labels = compute_qg10_labels(s)

        amap = compute_amap(s, labels)

        # Externals
        @test amap[1, 1] == 1
        @test amap[2, 1] == 2
        @test amap[3, 1] == 3
        @test amap[4, 1] == 4
        # Internal vertex 5: slots 1,2 from externals (back-write); slot 3 internal
        @test amap[5, 1] == 1
        @test amap[5, 2] == 2
        @test amap[5, 3] == 5      # the internal 5-6 edge
        # Internal vertex 6: slots 1,2 from externals; slot 3 internal
        @test amap[6, 1] == 3
        @test amap[6, 2] == 4
        @test amap[6, 3] == 5      # same id ⇒ pairing invariant
    end

    @testset "φ³ self-energy 1L: parallel edges get distinct ids" begin
        # 2 ext (1,2), 2 int (3,4) with TWO parallel 3-4 edges.
        # ege[3,4] = nleg+1 = 3. Parallel edges occupy ids 3 and 4 (=ege+0/+1).
        p = Partition(Int8(2), Int8[2], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 3] = 1
        s.xg[2, 4] = 1
        s.xg[3, 4] = 2     # parallel edges
        s.xn[3]    = 1
        s.xn[4]    = 1
        labels = compute_qg10_labels(s)

        amap = compute_amap(s, labels)

        @test amap[1, 1] == 1
        @test amap[2, 1] == 2
        @test amap[3, 1] == 1     # ext back-write
        @test amap[4, 1] == 2
        @test amap[3, 2] == 3     # first parallel id (= ege[3,4])
        @test amap[3, 3] == 4     # second parallel id
        @test amap[4, 2] == 3     # pairs with amap[3,2]
        @test amap[4, 3] == 4
    end

    @testset "φ³ tadpole 1L: self-loop both slots have same id" begin
        # 1 ext (1), 1 int (2) with self-loop.
        # ege[2,2] = nleg+1 = 2. Self-loop's two half-edge slots both get id 2.
        p = Partition(Int8(1), Int8[1], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 2] = 1
        s.xg[2, 2] = 2     # self-loop
        s.xn[2]    = 1
        labels = compute_qg10_labels(s)

        amap = compute_amap(s, labels)

        @test amap[1, 1] == 1
        @test amap[2, 1] == 1   # ext back-write
        @test amap[2, 2] == 2   # self-loop slot 1
        @test amap[2, 3] == 2   # self-loop slot 2 — SAME id (one edge)
    end

    @testset "Pairing invariant: amap[v,s] == amap[vmap[v,s], lmap[v,s]]" begin
        # Comprehensive invariant check across all the above topologies.
        # If this passes, downstream consumers can treat amap as edge id
        # without worrying about endpoint orientation.
        for (descr, build) in [
            ("φ→φφ tree", () -> begin
                p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
                s = TopoState(p)
                s.xg[1,4]=1; s.xg[2,4]=1; s.xg[3,4]=1
                s.xn[4]=3
                s
            end),
            ("ee→μμ s-channel", () -> begin
                p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
                s = TopoState(p)
                s.xg[1,5]=1; s.xg[2,5]=1; s.xg[3,6]=1; s.xg[4,6]=1; s.xg[5,6]=1
                s.xn[5]=2; s.xn[6]=2
                s
            end),
            ("φ³ 1L bubble", () -> begin
                p = Partition(Int8(2), Int8[2], Int8(3), Int8(1))
                s = TopoState(p)
                s.xg[1,3]=1; s.xg[2,4]=1; s.xg[3,4]=2
                s.xn[3]=1; s.xn[4]=1
                s
            end),
            ("φ³ tadpole", () -> begin
                p = Partition(Int8(1), Int8[1], Int8(3), Int8(1))
                s = TopoState(p)
                s.xg[1,2]=1; s.xg[2,2]=2
                s.xn[2]=1
                s
            end),
        ]
            @testset "$descr" begin
                s = build()
                labels = compute_qg10_labels(s)
                amap = compute_amap(s, labels)
                for v in 1:Int(s.n)
                    for slot in 1:Int(s.vdeg[v])
                        peer_v = Int(labels.vmap[v, slot])
                        peer_s = Int(labels.lmap[v, slot])
                        @test amap[v, slot] == amap[peer_v, peer_s]
                    end
                end
            end
        end
    end

end
