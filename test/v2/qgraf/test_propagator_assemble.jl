#  Phase 18a-3: per-edge propagator factors.
#
#  build_propagators(state, labels, pmap, edge_mom, model) returns
#  one Propagator per internal edge with field, momentum, num and denom
#  factors.  Tree-only — no PaVe machinery.
#
#  Numerator convention (matches existing src/v2/amplitude.jl):
#    Boson  : alg(1)  (metric is implicit via Lorentz-index sharing at vertices)
#    Fermion: DiracExpr(p-slash + m) via propagator_num
#    Scalar : alg(1)  via propagator_num
#  Denominator: AlgSum representing (p² - m²) = pair(mom, mom) - m².

using Test

using Feynfeld
using Feynfeld.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels, route_momenta, compute_amap,
                            build_propagators, Propagator

@testset "Phase 18a-3: per-edge propagator factors" begin

    @testset "ee→μμ s-channel: photon propagator" begin
        # 4 ext (1=e, 2=e_bar, 3=mu, 4=mu_bar), 2 int (5,6).
        # Edge 5-6 is the photon (γ exchange).
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

        # pmap: external slots + internal photon
        pmap = fill(:_, 6, MAX_V)
        pmap[1, 1] = :e;       pmap[5, 1] = :e_bar
        pmap[2, 1] = :e_bar;   pmap[5, 2] = :e
        pmap[3, 1] = :mu;      pmap[6, 1] = :mu_bar
        pmap[4, 1] = :mu_bar;  pmap[6, 2] = :mu
        pmap[5, 3] = :gamma;   pmap[6, 3] = :gamma

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edge_mom = route_momenta(s, labels, ext_moms)

        propagators = build_propagators(s, labels, pmap, edge_mom, qed_model())

        @test length(propagators) == 1
        p_γ = propagators[1]
        @test p_γ.field == :gamma
        @test (p_γ.v_lo, p_γ.v_hi, p_γ.parallel_idx) == (5, 6, 1)
        # Photon carries p1+p2 in qgraf "all incoming" convention
        @test p_γ.mom == momentum_sum([(1//1, Momentum(:p1)),
                                       (1//1, Momentum(:p2))])
        # Boson numerator: implicit metric, returned as scalar 1
        @test p_γ.num == alg(1)
        # Photon is massless ⇒ denominator = (p1+p2)² = pair(mom, mom)
        @test p_γ.denom == alg(pair(p_γ.mom, p_γ.mom))
    end

    @testset "φ³ φφ→φφ s-channel: scalar propagator" begin
        # 4 ext (1..4 = :phi), 2 int (5,6).  Internal edge is φ.
        # Like ee→μμ topologically; species drives propagator branch.
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

        propagators = build_propagators(s, labels, pmap, edge_mom,
                                         phi3_model(mass=:zero))

        @test length(propagators) == 1
        p_φ = propagators[1]
        @test p_φ.field == :phi
        @test p_φ.num == alg(1)   # scalar numerator = 1 (propagator_num)
        # φ³ φ is massless in phi3_model → denom = pair(p1+p2, p1+p2)
        expected_mom = momentum_sum([(1//1, Momentum(:p1)),
                                     (1//1, Momentum(:p2))])
        @test p_φ.mom   == expected_mom
        @test p_φ.denom == alg(pair(expected_mom, expected_mom))
    end

    @testset "φ³ self-energy 1L bubble: two parallel scalar propagators" begin
        # 2 ext (1,2 = :phi), 2 int (3,4) with TWO parallel 3-4 edges.
        # Expected: 2 propagators, both :phi, distinct momenta (p1-k and k).
        p = Partition(Int8(2), Int8[2], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 3] = 1
        s.xg[2, 4] = 1
        s.xg[3, 4] = 2
        s.xn[3]    = 1
        s.xn[4]    = 1
        labels = compute_qg10_labels(s)

        pmap = fill(:_, 4, MAX_V)
        for v in 1:4, slot in 1:Int(s.vdeg[v])
            pmap[v, slot] = :phi
        end

        ext_moms  = [Momentum(:p1), Momentum(:p2)]
        loop_moms = [Momentum(:k)]
        edge_mom  = route_momenta(s, labels, ext_moms; loop_moms)

        propagators = build_propagators(s, labels, pmap, edge_mom,
                                         phi3_model(mass=:zero))

        @test length(propagators) == 2
        @test all(p -> p.field == :phi, propagators)
        @test all(p -> p.num == alg(1), propagators)
        @test propagators[1].parallel_idx == 1
        @test propagators[2].parallel_idx == 2
        # Distinct momenta from the leaf-peel result
        moms = [p.mom for p in propagators]
        @test momentum_sum([(1//1, Momentum(:p1)),
                            (-1//1, Momentum(:k))]) in moms
        @test Momentum(:k) in moms
    end

    @testset "φ³ tadpole 1L: self-loop scalar propagator" begin
        # 1 ext (1), 1 int (2) with self-loop. Expected: 1 propagator,
        # field=:phi, momentum = loop momentum k.
        p = Partition(Int8(1), Int8[1], Int8(3), Int8(1))
        s = TopoState(p)
        s.xg[1, 2] = 1
        s.xg[2, 2] = 2
        s.xn[2]    = 1
        labels = compute_qg10_labels(s)

        pmap = fill(:_, 2, MAX_V)
        for v in 1:2, slot in 1:Int(s.vdeg[v])
            pmap[v, slot] = :phi
        end

        ext_moms = [Momentum(:p1)]
        edge_mom = route_momenta(s, labels, ext_moms)   # default :k1

        propagators = build_propagators(s, labels, pmap, edge_mom,
                                         phi3_model(mass=:zero))

        @test length(propagators) == 1
        p_self = propagators[1]
        @test p_self.field == :phi
        @test (p_self.v_lo, p_self.v_hi) == (2, 2)
        @test p_self.mom   == Momentum(:k1)
        @test p_self.num   == alg(1)
        @test p_self.denom == alg(pair(Momentum(:k1), Momentum(:k1)))
    end

    # Phase 18b-2 (feynfeld-h3pb): fermion propagator numerator must accept
    # composite momentum. γ^μ is linear in its argument, so for p = Σ cᵢ pᵢ
    # the numerator (p̸ + m) expands to Σ cᵢ p̸ᵢ + m·I.
    @testset "fermion propagator numerator with composite momentum" begin
        _propnum = Feynfeld.QgrafPort._propagator_numerator
        p1 = Momentum(:p1)
        p2 = Momentum(:p2)

        @testset "massless: p̸1 + p̸2" begin
            mom = momentum_sum([(1//1, p1), (1//1, p2)])
            @test mom isa MomentumSum
            expected = DiracExpr(DiracChain([GS(p1)])) +
                       DiracExpr(DiracChain([GS(p2)]))
            @test _propnum(Fermion(), mom, 0//1) == expected
        end

        @testset "massive: p̸1 + p̸2 + m·I" begin
            mom = momentum_sum([(1//1, p1), (1//1, p2)])
            mass = 1//1
            expected = DiracExpr(DiracChain([GS(p1)])) +
                       DiracExpr(DiracChain([GS(p2)])) +
                       mass * DiracExpr(alg(1))
            @test _propnum(Fermion(), mom, mass) == expected
        end

        @testset "negative coefficient: p̸1 − p̸2" begin
            mom = momentum_sum([(1//1, p1), (-1//1, p2)])
            expected = DiracExpr(DiracChain([GS(p1)])) -
                       DiracExpr(DiracChain([GS(p2)]))
            @test _propnum(Fermion(), mom, 0//1) == expected
        end

        @testset "rational coefficient: ½ p̸1" begin
            # Single-term momentum_sum with non-unit coefficient stays a
            # MomentumSum (won't degenerate to bare Momentum).
            mom = momentum_sum([(1//2, p1)])
            @test mom isa MomentumSum
            expected = (1//2) * DiracExpr(DiracChain([GS(p1)]))
            @test _propnum(Fermion(), mom, 0//1) == expected
        end
    end

    # Phase 18b-2 (feynfeld-h3pb): build_propagators must succeed on a
    # Compton-like topology where the internal edge is an off-shell fermion
    # carrying p1+p2.
    @testset "Compton-like s-channel: internal fermion propagator" begin
        # Topology: 4 ext (1=e, 2=γ, 3=e, 4=γ), 2 int (5,6). Edge 5-6 is
        # the off-shell electron. Mirrors the photon-exchange topology
        # above but with the internal field flipped to an electron.
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

        ext_moms = [Momentum(:p1), Momentum(:p2), Momentum(:p3), Momentum(:p4)]
        edge_mom = route_momenta(s, labels, ext_moms)

        propagators = build_propagators(s, labels, pmap, edge_mom, qed_model())

        @test length(propagators) == 1
        p_e = propagators[1]
        @test p_e.field == :e
        @test (p_e.v_lo, p_e.v_hi, p_e.parallel_idx) == (5, 6, 1)
        # All-incoming convention: internal fermion carries p1+p2
        expected_mom = momentum_sum([(1//1, Momentum(:p1)),
                                     (1//1, Momentum(:p2))])
        @test p_e.mom == expected_mom
        # Electron has placeholder mass 1//1 in QED model (m_e ≠ :zero).
        # Numerator: p̸1 + p̸2 + 1·I.
        expected_num = DiracExpr(DiracChain([GS(Momentum(:p1))])) +
                       DiracExpr(DiracChain([GS(Momentum(:p2))])) +
                       (1//1) * DiracExpr(alg(1))
        @test p_e.num == expected_num
        # Denominator: (p1+p2)² − m² (placeholder mass 1//1).
        @test p_e.denom == alg(pair(expected_mom, expected_mom)) - (1//1) * alg(1)
    end

end
