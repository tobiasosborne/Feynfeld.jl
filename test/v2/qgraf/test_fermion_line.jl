#  Phase 18a-6: fermion-line traversal.
#
#  walk_fermion_lines(state, labels, pmap, physical_moms, n_inco, model)
#  returns one FermionLine per fermion path through the diagram. For
#  tree QED 2→2, each fermion line touches exactly ONE internal vertex
#  (no internal fermion propagator); 1-loop and Compton cases (multi-
#  vertex fermion lines) are deferred to Phase 18b.

using Test

using Feynfeld
using Feynfeld.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels,
                            walk_fermion_lines, FermionLine

@testset "Phase 18a-6: fermion-line traversal" begin

    @testset "ee→μμ s-channel: 2 lines, one per vertex" begin
        # Vertex 5 hosts the e-line (vbar(p2) γ^μ u(p1));
        # vertex 6 hosts the μ-line (ubar(k1) γ_μ v(k2)).
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

        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]
        lines = walk_fermion_lines(s, labels, pmap, physical_moms, 2,
                                    qed_model(m_e=:zero, m_mu=:zero))

        @test length(lines) == 2

        e_line  = lines[findfirst(l -> l.vertex == 5, lines)]
        mu_line = lines[findfirst(l -> l.vertex == 6, lines)]

        # e-line: ext leg 2 = e_bar (incoming antiparticle) → vbar at bar end;
        #         ext leg 1 = e (incoming particle)         → u    at plain end.
        @test e_line.bar_leg   == 2
        @test e_line.plain_leg == 1
        # bar slot at vertex 5 connects via vmap to ext leg 2 (the bar leg)
        @test Int(labels.vmap[e_line.vertex, e_line.bar_slot])   == 2
        @test Int(labels.vmap[e_line.vertex, e_line.plain_slot]) == 1

        # μ-line: ext leg 3 = mu (outgoing particle)     → ubar at bar end;
        #         ext leg 4 = mu_bar (outgoing antiparticle) → v at plain end.
        @test mu_line.bar_leg   == 3
        @test mu_line.plain_leg == 4
        @test Int(labels.vmap[mu_line.vertex, mu_line.bar_slot])   == 3
        @test Int(labels.vmap[mu_line.vertex, mu_line.plain_slot]) == 4
    end

    @testset "φ³ φφ→φφ: no fermion lines" begin
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

        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]
        lines = walk_fermion_lines(s, labels, pmap, physical_moms, 2,
                                    phi3_model(mass=:zero))
        @test isempty(lines)
    end

end
