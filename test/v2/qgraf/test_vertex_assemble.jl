#  Phase 18a-4: per-vertex Lorentz factors.
#
#  build_vertices(state, labels, pmap, edge_mom, model) -> Dict{Int, DiracExpr}
#  emits one vertex factor per internal vertex, allocating a shared
#  Lorentz index per attached boson edge so factors at the two endpoints
#  contract automatically when chains multiply.
#
#  Convention (matches src/v2/amplitude.jl): the Lorentz index for an
#  edge is `LorentzIndex(Symbol(:mu_l_, edge_id), DimD())` where edge_id
#  is the global half-edge id from compute_amap.

using Test

include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels, route_momenta,
                            build_vertices

@testset "Phase 18a-4: per-vertex Lorentz factors" begin

    @testset "ee→μμ s-channel: γ^μ at both vertices, shared index" begin
        # 4 ext, 2 int (5,6).  Both internal vertices are QED 3-vertices
        # with a photon attached at slot 3.  The photon edge has amap=5
        # (= rhop1 = nleg+1) → both vertices use index :mu_l_5.
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

        pmap = fill(:_, 6, MAX_V)
        pmap[1, 1] = :e;       pmap[5, 1] = :e_bar
        pmap[2, 1] = :e_bar;   pmap[5, 2] = :e
        pmap[3, 1] = :mu;      pmap[6, 1] = :mu_bar
        pmap[4, 1] = :mu_bar;  pmap[6, 2] = :mu
        pmap[5, 3] = :gamma;   pmap[6, 3] = :gamma

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edge_mom = route_momenta(s, labels, ext_moms)

        vertices = build_vertices(s, labels, pmap, edge_mom, qed_model())

        @test length(vertices) == 2
        @test haskey(vertices, 5) && haskey(vertices, 6)

        # Both vertices are γ^{:mu_l_5} (same shared index ⇒ contraction
        # is automatic when chains multiply at amplitude assembly time).
        mu = LorentzIndex(Symbol(:mu_l_, 5), DimD())
        expected = DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))
        @test vertices[5] == expected
        @test vertices[6] == expected
    end

    @testset "φ³ φφ→φφ s-channel: scalar vertices (no Lorentz structure)" begin
        # Same shape, but all fields are :phi → no Lorentz index, no Dirac.
        # Vertex factor at each internal vertex is just DiracExpr(alg(1)).
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

        pmap = fill(:_, 6, MAX_V)
        for v in 1:6, slot in 1:Int(s.vdeg[v])
            pmap[v, slot] = :phi
        end

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edge_mom = route_momenta(s, labels, ext_moms)

        vertices = build_vertices(s, labels, pmap, edge_mom, phi3_model())
        @test length(vertices) == 2
        @test vertices[5] == DiracExpr(alg(1))
        @test vertices[6] == DiracExpr(alg(1))
    end

end
