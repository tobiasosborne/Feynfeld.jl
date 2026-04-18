#  Phase 18a-1: leaf-peel momentum routing.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13400-13559
#  Cross-ref: refs/qgraf/ALGORITHM.md §5.2
#
#  qgraf convention: every external momentum flows INTO the diagram. Caller
#  negates outgoing-physical momenta if physical-flow sign is desired.

using Test

using Feynfeld: Momentum, MomentumSum, momentum_sum
using Feynfeld.QgrafPort: Partition, TopoState, qg21_enumerate!,
                            compute_qg10_labels, route_momenta, InternalEdge

@testset "Phase 18a-1: leaf-peel momentum routing" begin

    @testset "φ→φφ tree (1 internal vertex, 0 internal edges)" begin
        # 3 externals (vertices 1,2,3) all attach to one internal vertex (4).
        # No internal edges → leaf-peel produces empty result.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 4] = 1
        s.xg[2, 4] = 1
        s.xg[3, 4] = 1
        s.xn[4]    = 3
        labels = compute_qg10_labels(s)

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3)]
        edges = route_momenta(s, labels, ext_moms)

        @test edges.n_ext == 3
        @test edges.ext_moms == ext_moms
        @test isempty(edges.internal)
    end

    @testset "ee→μμ t-channel: internal edge carries p1+p3" begin
        # Same shape as the s-channel test, but externals attach differently:
        # 1,3 attach at v=5; 2,4 attach at v=6.  Expected (qgraf
        # "all incoming"): edge carries p1+p3.  (Physical convention p1-k1
        # is recovered if the caller passes Momentum(:k1) negated as p3.)
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 5] = 1
        s.xg[3, 5] = 1
        s.xg[2, 6] = 1
        s.xg[4, 6] = 1
        s.xg[5, 6] = 1
        s.xn[5]    = 2
        s.xn[6]    = 2
        labels = compute_qg10_labels(s)

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edges = route_momenta(s, labels, ext_moms)

        @test length(edges.internal) == 1
        e = edges.internal[1]
        @test e.momentum == momentum_sum([(1//1, Momentum(:p1)),
                                          (1//1, Momentum(:p3))])
        @test e.edge_type == :rb
    end

    @testset "ee→μμ s-channel: internal edge carries p1+p2" begin
        # 4 ext (1..4), 2 int (5,6).  Externals 1,2 attach at v=5 (incoming
        # side); externals 3,4 attach at v=6 (outgoing side).  Single
        # internal edge 5-6 in the spanning tree, no chords (tree topology).
        # Expected (qgraf "all incoming" convention): edge (5,6) carries
        # p1 + p2 (sum of externals at the lower-indexed endpoint).
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

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edges = route_momenta(s, labels, ext_moms)

        @test length(edges.internal) == 1
        e = edges.internal[1]
        @test (e.v_lo, e.v_hi, e.parallel_idx) == (5, 6, 1)
        @test e.momentum == momentum_sum([(1//1, Momentum(:p1)),
                                          (1//1, Momentum(:p2))])
        @test e.edge_type == :rb   # regular bridge: only external momenta
    end

    @testset "φ³ self-energy 1L: bubble (1 chord, head-match sign flip)" begin
        # 2 ext (1,2), 2 int (3,4).  Single edge 1-3, single edge 2-4,
        # TWO parallel edges 3-4 (bubble).  Spanning tree picks the first
        # 3-4 copy; the second is a chord carrying loop momentum k.
        # Expected (qgraf "all incoming"):
        #   internal edge (3,4,1) tree branch: p1 - k    (head match → subtract)
        #   internal edge (3,4,2) chord       : k
        p = Partition(Int8(2), Int8[2], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 3] = 1
        s.xg[2, 4] = 1
        s.xg[3, 4] = 2
        s.xn[3]    = 1
        s.xn[4]    = 1
        labels = compute_qg10_labels(s)

        ext_moms  = [Momentum(:p1), Momentum(:p2)]
        loop_moms = [Momentum(:k)]
        edges = route_momenta(s, labels, ext_moms; loop_moms)

        @test length(edges.internal) == 2
        e_tree, e_chord = edges.internal[1], edges.internal[2]
        @test (e_tree.v_lo, e_tree.v_hi, e_tree.parallel_idx)  == (3, 4, 1)
        @test (e_chord.v_lo, e_chord.v_hi, e_chord.parallel_idx) == (3, 4, 2)

        @test e_tree.momentum  == momentum_sum([(1//1, Momentum(:p1)),
                                                (-1//1, Momentum(:k))])
        @test e_chord.momentum == Momentum(:k)

        @test e_tree.edge_type  == :rnb   # carries loop momentum, not a self-loop
        @test e_chord.edge_type == :rnb
    end

    @testset "φ³ tadpole 1L: self-loop chord (snb edge type)" begin
        # 1 ext (1), 1 int (2) with a self-loop.  vdeg[2] = 1 (ext) + 2
        # (self-loop half-edges) = 3.  No internal tree edges, so leaf-peel
        # never fires; the only internal edge is the self-loop chord.
        p = Partition(Int8(1), Int8[1], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 2] = 1
        s.xg[2, 2] = 2   # qgraf stores self-loops as 2 × #self-loops
        s.xn[2]    = 1
        labels = compute_qg10_labels(s)

        ext_moms = [Momentum(:p1)]
        edges = route_momenta(s, labels, ext_moms)   # default loop name :k1

        @test length(edges.internal) == 1
        e = edges.internal[1]
        @test (e.v_lo, e.v_hi, e.parallel_idx) == (2, 2, 1)
        @test e.momentum  == Momentum(:k1)
        @test e.edge_type == :snb   # self-loop carrying loop momentum
    end

    @testset "qg21_enumerate integration: 4-ext deg-3 tree topology" begin
        # qg21 emits ONE topology for partition (4, [2], 3, 0) up to iso —
        # the s/t/u channels arise from field assignment, not topology.
        # Verify route_momenta runs end-to-end on a state qg21 actually
        # produces (not a hand-built one).
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        @test length(emits) == 1   # canonical topology unique
        state = emits[1]
        labels = compute_qg10_labels(state)
        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edges = route_momenta(state, labels, ext_moms)
        # 1 internal edge between the two deg-3 internal vertices
        @test length(edges.internal) == 1
        e = edges.internal[1]
        @test e.edge_type == :rb
        # Internal edge should carry sum of two consecutive externals (qg21
        # canonicalisation puts {1,2} at one int vertex, {3,4} at the other).
        @test e.momentum == momentum_sum([(1//1, Momentum(:p1)),
                                          (1//1, Momentum(:p2))])
    end

end
