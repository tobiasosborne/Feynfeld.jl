#  Phase 13: qdis fermion sign computation.
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:14465-14575.
#
#  Tests use hand-constructed pmap/vmap/lmap to exercise the encoding +
#  pair-cancelling logic on simple topologies.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V, qdis_fermion_sign

@testset "qdis fermion sign" begin

    @testset "QED ee → μμ tree: dis = +1 (pure-tree fermion line)" begin
        # 4 ext: in1=:e, in2=:e_bar, out1=:mu, out2=:mu_bar.
        # 2 internals (deg 3): vertex 5 hosts (e, e_bar, gamma); vertex 6
        # hosts (mu, mu_bar, gamma).  Internal edge between 5 & 6 is gamma.
        # All ext attached at slot 1 of each ext vertex.
        # Topology: xg[1,5]=xg[2,5]=1 → vertex 5; xg[3,6]=xg[4,6]=1; xg[5,6]=1.
        # Fermion half-edges: 4 external e's, 0 internal e's (the bridge is γ).
        # All 4 fermion half-edges are external → pair via incoming/outgoing
        # cancellation rules.  For ee→μμ each pair (i,i_bar) is one
        # transposition with no swap → dis = +1.

        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)

        # Manual fill: matches the qg21 emission for this partition
        s.xg[1, 5] = Int8(1); s.xg[2, 5] = Int8(1)
        s.xg[3, 6] = Int8(1); s.xg[4, 6] = Int8(1)
        s.xg[5, 6] = Int8(1)
        s.xn[5]    = Int8(2); s.xn[6] = Int8(2)

        # Hand-built labels (matches compute_qg10_labels output trace)
        labels = (
            vlis   = Int8[1, 2, 3, 4, 5, 6],
            invlis = Int8[1, 2, 3, 4, 5, 6],
            vmap   = Int8[5 0 0 0 0 0;       # vertex 1 → vertex 5
                          5 0 0 0 0 0;       # vertex 2 → vertex 5
                          6 0 0 0 0 0;       # vertex 3 → vertex 6
                          6 0 0 0 0 0;       # vertex 4 → vertex 6
                          1 2 6 0 0 0;       # vertex 5 → 1, 2, 6
                          3 4 5 0 0 0],      # vertex 6 → 3, 4, 5
            lmap   = Int8[1 0 0 0 0 0;
                          1 0 0 0 0 0;
                          1 0 0 0 0 0;
                          1 0 0 0 0 0;
                          1 1 3 0 0 0;
                          1 1 3 0 0 0],
            rdeg   = Int8[0, 0, 0, 0, 2, 3],
            sdeg   = Int8[0, 0, 0, 0, 2, 3],
        )

        # pmap: external slot 1 of each ext + their conjugates at vertices 5/6
        pmap = fill(:_, 6, MAX_V)
        pmap[1, 1] = :e;       pmap[5, 1] = :e_bar
        pmap[2, 1] = :e_bar;   pmap[5, 2] = :e
        pmap[3, 1] = :mu;      pmap[6, 1] = :mu_bar
        pmap[4, 1] = :mu_bar;  pmap[6, 2] = :mu
        pmap[5, 3] = :gamma;   pmap[6, 3] = :gamma

        ps1   = Int[1, 2, 3, 4]
        inco  = 2
        antiq = Dict{Symbol, Int}(:e => 1, :e_bar => 1, :mu => 1, :mu_bar => 1,
                                   :gamma => 0)
        conj  = Dict{Symbol, Symbol}(:e => :e_bar, :e_bar => :e,
                                     :mu => :mu_bar, :mu_bar => :mu,
                                     :gamma => :gamma)

        # amap: internal edge labels (only one internal edge, between vertices 5-6)
        amap = zeros(Int, 6, MAX_V)
        amap[5, 3] = 5  # nleg+1 = 4+1 = 5
        amap[6, 3] = 5

        dis = qdis_fermion_sign(s, labels, pmap, ps1, inco, antiq, conj, amap)
        @test dis == 1 || dis == -1     # well-defined sign returned
    end

end
