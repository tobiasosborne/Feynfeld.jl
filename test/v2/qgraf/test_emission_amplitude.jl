#  Phase 18a-7: master assembler emission_to_amplitude.
#
#  Composes the per-piece builders (route_momenta, compute_amap,
#  build_propagators, build_vertices, build_externals, walk_fermion_lines)
#  into an AmplitudeBundle for one emission.  Tree-only.

using Test

using Feynfeld
using Feynfeld.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels,
                            emission_to_amplitude, AmplitudeBundle

@testset "Phase 18a-7: emission_to_amplitude" begin

    @testset "ee→μμ s-channel: bundle has 2-line amplitude + photon denom" begin
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

        ps1 = Int[1, 2, 3, 4]
        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]

        bundle = emission_to_amplitude(s, labels, ps1, pmap,
                                        qed_model(m_e=:zero, m_mu=:zero);
                                        physical_moms, n_inco=2)

        @test bundle isa AmplitudeBundle
        @test bundle.amplitude isa DiracExpr

        # Photon denom (p1+p2)² since the two ext_moms at the s-channel
        # vertex are p1 and p2.  qgraf "all incoming": ext_moms[3,4] are
        # negated outgoing — so route_momenta uses [p1, p2, -k1, -k2]
        # internally; the photon edge carries p1 + p2.
        s_channel_mom = momentum_sum([(1//1, Momentum(:p1)),
                                      (1//1, Momentum(:p2))])
        @test bundle.denoms == [alg(pair(s_channel_mom, s_channel_mom))]

        # qdis convention is well-defined ±1 (pinned at validation time
        # in Phase 18a-9 against the hand-built path).
        @test bundle.fermion_sign in (-1, 1)

        # Tree-level QED: no automorphisms → S_local = 1.
        @test bundle.sym_factor == 1//1
    end

    # Phase 18b-3 (feynfeld-a7f2): Compton tree s-channel exercises a
    # multi-vertex fermion line with one internal fermion propagator.
    @testset "Compton tree s-channel: bundle amplitude ≡ handbuilt" begin
        # 4 ext (1=e_in, 2=γ_in, 3=e_out, 4=γ_out), 2 int (5,6).
        # Internal edge 5-6 = off-shell electron carrying p1+k1.
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
        pmap[2, 1] = :gamma;   pmap[5, 2] = :gamma
        pmap[3, 1] = :e_bar;   pmap[6, 1] = :e
        pmap[4, 1] = :gamma;   pmap[6, 2] = :gamma
        pmap[5, 3] = :e;       pmap[6, 3] = :e_bar

        ps1 = Int[1, 2, 3, 4]
        physical_moms = [Momentum(:p1), Momentum(:k1),
                         Momentum(:p2), Momentum(:k2)]
        # Both electrons are physical particles; phys_anti needed because
        # qgraf labels the outgoing electron as :e_bar (all-incoming flip).
        phys_anti = [false, false, false, false]

        bundle = emission_to_amplitude(s, labels, ps1, pmap,
                                        qed_model(m_e=:zero, m_mu=:zero);
                                        physical_moms, n_inco=2,
                                        phys_anti=phys_anti)

        @test bundle isa AmplitudeBundle
        @test length(bundle.line_chains) == 1
        # One internal propagator → denoms has length 1 = (p1+k1)².
        s_channel_mom = momentum_sum([(1//1, Momentum(:p1)),
                                      (1//1, Momentum(:k1))])
        @test bundle.denoms == [alg(pair(s_channel_mom, s_channel_mom))]

        # Hand-built reference using the qgraf-port :mu_l_<edge_id>
        # naming convention. The line walks bar(leg 3)→v6→prop→v5→plain(leg 1).
        # mu_l_4 at v6 (boson edge to ext 4); mu_l_2 at v5 (boson edge to ext 2).
        # Propagator num: (p̸1 + k̸1) since m_e=:zero.
        mu_l_2 = LorentzIndex(:mu_l_2, DimD())
        mu_l_4 = LorentzIndex(:mu_l_4, DimD())
        ubar_p2 = ubar(Momentum(:p2))
        u_p1    = u(Momentum(:p1))
        gamma2  = DiracExpr(DiracChain([DiracGamma(LISlot(mu_l_2))]))
        gamma4  = DiracExpr(DiracChain([DiracGamma(LISlot(mu_l_4))]))
        prop    = DiracExpr(DiracChain([GS(Momentum(:p1))])) +
                  DiracExpr(DiracChain([GS(Momentum(:k1))]))
        expected = DiracExpr(Feynfeld.dot(ubar_p2)) * gamma4 * prop * gamma2 *
                   DiracExpr(Feynfeld.dot(u_p1))

        @test bundle.amplitude == expected
        @test bundle.line_chains[1] == expected
    end

    @testset "φ³ φφ→φφ s-channel: scalar bundle (no fermion lines)" begin
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

        ps1 = Int[1, 2, 3, 4]
        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]

        bundle = emission_to_amplitude(s, labels, ps1, pmap,
                                        phi3_model(mass=:zero);
                                        physical_moms, n_inco=2)

        @test bundle.amplitude == DiracExpr(alg(1))   # no Dirac structure
        s_channel_mom = momentum_sum([(1//1, Momentum(:p1)),
                                      (1//1, Momentum(:p2))])
        @test bundle.denoms == [alg(pair(s_channel_mom, s_channel_mom))]
        @test bundle.fermion_sign == 1                # no fermions
        @test bundle.sym_factor == 1//1
    end

end
