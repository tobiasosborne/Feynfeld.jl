#  Phase 18b-3: multi-vertex fermion-line traversal.
#
#  walk_fermion_lines(state, labels, pmap, physical_moms, n_inco, model)
#  returns one FermionLine per fermion path through the diagram. For
#  tree QED 2→2 each line has a single internal vertex (no propagator);
#  Compton tree lines have two vertices joined by one internal fermion
#  propagator (Phase 18b-3 = feynfeld-a7f2).

using Test

using Feynfeld
using Feynfeld.QgrafPort: Partition, TopoState, MAX_V,
                            compute_qg10_labels,
                            walk_fermion_lines, FermionLine

@testset "Phase 18b-3: fermion-line traversal" begin

    @testset "ee→μμ s-channel: 2 single-vertex lines" begin
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

        e_line  = lines[findfirst(l -> l.vertices == [5], lines)]
        mu_line = lines[findfirst(l -> l.vertices == [6], lines)]

        # Single-vertex tree lines: no internal fermion propagators.
        @test isempty(e_line.propagator_edge_ids)
        @test isempty(mu_line.propagator_edge_ids)

        # e-line: ext leg 2 = e_bar (incoming antiparticle) → vbar at bar end;
        #         ext leg 1 = e (incoming particle)         → u    at plain end.
        @test e_line.bar_leg   == 2
        @test e_line.plain_leg == 1
        # bar slot at vertex 5 connects via vmap to ext leg 2 (the bar leg)
        @test Int(labels.vmap[e_line.vertices[1], e_line.in_slots[1]])  == 2
        @test Int(labels.vmap[e_line.vertices[1], e_line.out_slots[1]]) == 1

        # μ-line: ext leg 3 = mu (outgoing particle)     → ubar at bar end;
        #         ext leg 4 = mu_bar (outgoing antiparticle) → v at plain end.
        @test mu_line.bar_leg   == 3
        @test mu_line.plain_leg == 4
        @test Int(labels.vmap[mu_line.vertices[1], mu_line.in_slots[1]])  == 3
        @test Int(labels.vmap[mu_line.vertices[1], mu_line.out_slots[1]]) == 4
    end

    @testset "Compton tree s-channel: 1 line, 2 vertices, 1 propagator" begin
        # Physical: e(p1) + γ(k1) → e(p2) + γ(k2).
        # Topology: 4 ext (1=e_in, 2=γ_in, 3=e_out, 4=γ_out), 2 int (5,6).
        # Internal edge 5-6 is the off-shell electron carrying p1+k1.
        # Vertex 5 hosts ext legs 1,2; vertex 6 hosts ext legs 3,4.
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
        # qgraf "all incoming": physical incoming e → :e, incoming γ → :gamma,
        # physical outgoing e → :e_bar (CP flip), outgoing γ → :gamma.
        pmap[1, 1] = :e;       pmap[5, 1] = :e_bar
        pmap[2, 1] = :gamma;   pmap[5, 2] = :gamma
        pmap[3, 1] = :e_bar;   pmap[6, 1] = :e
        pmap[4, 1] = :gamma;   pmap[6, 2] = :gamma
        pmap[5, 3] = :e;       pmap[6, 3] = :e_bar

        physical_moms = [Momentum(:p1), Momentum(:k1),
                         Momentum(:p2), Momentum(:k2)]
        # Both electrons are physical particles (not antiparticles); the
        # qgraf :e_bar label on outgoing leg 3 reflects the all-incoming
        # convention, not physical anti-status. Pass phys_anti so that
        # _spinor_dispatch sees the right particle/antiparticle identity.
        phys_anti = [false, false, false, false]
        lines = walk_fermion_lines(s, labels, pmap, physical_moms, 2,
                                    qed_model(m_e=:zero, m_mu=:zero);
                                    phys_anti=phys_anti)

        @test length(lines) == 1
        line = lines[1]
        @test length(line.vertices) == 2
        @test length(line.in_slots) == 2
        @test length(line.out_slots) == 2
        @test length(line.propagator_edge_ids) == 1

        # The line walks v=5 → propagator → v=6 (or 6 → 5 depending on
        # which leg is :left). Either way the two internal vertices both
        # appear, and bar/plain legs are externals 1 and 3 (the e legs).
        @test Set(line.vertices) == Set([5, 6])
        @test Set([line.bar_leg, line.plain_leg]) == Set([1, 3])

        # Each step's in_slot must be a fermion slot that vmaps to either
        # the previous vertex/leg or the bar external, and similarly for
        # out_slot. Specifically, the propagator edge_id at the seam must
        # equal amap[v_k, out_slot_k] = amap[v_{k+1}, in_slot_{k+1}].
        amap = Feynfeld.QgrafPort.compute_amap(s, labels)
        v1, v2 = line.vertices[1], line.vertices[2]
        os1, is2 = line.out_slots[1], line.in_slots[2]
        @test amap[v1, os1] == line.propagator_edge_ids[1]
        @test amap[v2, is2] == line.propagator_edge_ids[1]
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
