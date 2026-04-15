#  Phase 18a-5: per-external spinor / polarisation factors.
#
#  build_externals(state, pmap, physical_moms, n_inco, model) returns
#  one ExternalFactor per external leg, holding (field, momentum,
#  in/out + antiparticle flags, spinor, position).  Mirrors the
#  src/v2/amplitude.jl `_spinor_and_position` dispatch.
#
#  Convention: physical_moms are PHYSICAL (not qgraf "all incoming")
#  — incoming legs and outgoing legs both pass their physical momenta.

using Test

include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                            build_externals, ExternalFactor

@testset "Phase 18a-5: per-external spinor/polarisation factors" begin

    @testset "ee→μμ massless: u, vbar, ubar, v dispatch" begin
        # 4 ext, with n_inco = 2.
        # leg 1: e   incoming particle      → u(p1)   right
        # leg 2: e_bar incoming antiparticle → vbar(p2) left
        # leg 3: mu  outgoing particle      → ubar(k1) left
        # leg 4: mu_bar outgoing antiparticle → v(k2)   right
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 5] = 1
        s.xg[2, 5] = 1
        s.xg[3, 6] = 1
        s.xg[4, 6] = 1
        s.xg[5, 6] = 1
        s.xn[5]    = 2
        s.xn[6]    = 2

        pmap = fill(:_, 6, MAX_V)
        pmap[1, 1] = :e
        pmap[2, 1] = :e_bar
        pmap[3, 1] = :mu
        pmap[4, 1] = :mu_bar
        # internal slots irrelevant for build_externals

        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]
        externals = build_externals(s, pmap, physical_moms, 2,
                                     qed_model(m_e=:zero, m_mu=:zero))

        @test length(externals) == 4

        e1, e2, e3, e4 = externals
        @test e1.field == :e       && e1.spinor == u(Momentum(:p1))    && e1.position == :right
        @test e2.field == :e_bar   && e2.spinor == vbar(Momentum(:p2)) && e2.position == :left
        @test e3.field == :mu      && e3.spinor == ubar(Momentum(:k1)) && e3.position == :left
        @test e4.field == :mu_bar  && e4.spinor == v(Momentum(:k2))    && e4.position == :right
    end

    @testset "φ³ φφ→φφ: scalar externals have no spinor / position" begin
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        s.xg[1, 5] = 1
        s.xg[2, 5] = 1
        s.xg[3, 6] = 1
        s.xg[4, 6] = 1
        s.xg[5, 6] = 1
        s.xn[5]    = 2
        s.xn[6]    = 2

        pmap = fill(:_, 6, MAX_V)
        for v in 1:6, slot in 1:Int(s.vdeg[v])
            pmap[v, slot] = :phi
        end

        physical_moms = [Momentum(:p1), Momentum(:p2),
                         Momentum(:k1), Momentum(:k2)]
        externals = build_externals(s, pmap, physical_moms, 2,
                                     phi3_model(mass=:zero))

        @test length(externals) == 4
        for e in externals
            @test e.field    == :phi
            @test e.spinor   === nothing
            @test e.position === nothing
        end
        # in/out flags still tracked (needed downstream for momentum routing)
        @test [e.incoming for e in externals] == [true, true, false, false]
    end

end
