#  Phase 17 audition: 3 dedup approaches must give the same qgraf-matching count.
#
#  Test cases (all hand-validated against qgraf golden masters):
#    phi3 tree φ→φφ          → 1
#    phi3 tree φφ→φφ          → 3   (s/t/u channels)
#    phi3 1L φ→φφ             → 7
#    phi3 1L φφ→φφ            → 39
#    QED ee→μμ tree            → 1
#    QED ee→μμ 1L (2-gen)      → 18
#
#  Each of (A) Burnside, (B) canonical-pmap, (C) pre-filter must agree.

using Test
using Feynfeld
using Feynfeld.QgrafPort: count_dedup_burnside, count_dedup_canonical,
                              count_dedup_prefilter

@testset "Phase 17 audition: dedup approaches agree" begin

    function _check(model, in_fields, out_fields; loops, expected)
        a = count_dedup_burnside(model,  in_fields, out_fields; loops=loops)
        b = count_dedup_canonical(model, in_fields, out_fields; loops=loops)
        c = count_dedup_prefilter(model, in_fields, out_fields; loops=loops)
        @test a == expected
        @test b == expected
        @test c == expected
    end

    @testset "phi3 tree φ→φφ → 1" begin
        _check(phi3_model(), [:phi], [:phi, :phi]; loops=0, expected=1)
    end

    @testset "phi3 tree φφ→φφ → 3 (s/t/u)" begin
        _check(phi3_model(), [:phi, :phi], [:phi, :phi]; loops=0, expected=3)
    end

    @testset "phi3 1L φ→φφ → 7" begin
        _check(phi3_model(), [:phi], [:phi, :phi]; loops=1, expected=7)
    end

    @testset "QED ee→μμ tree → 1" begin
        _check(qed_model(), [:e, :e], [:mu, :mu]; loops=0, expected=1)
    end

    @testset "phi3 1L φφ→φφ → 39 (qgraf golden master)" begin
        _check(phi3_model(), [:phi, :phi], [:phi, :phi]; loops=1, expected=39)
    end

    @testset "QED ee→μμ 1L (2-gen)" begin
        # Phase 12d (qgen.jl): _qgen_recurse enumerates all positional perms
        # of remaining slots + self-loop / multi-edge filters
        # (qgen:13921-13954). Session 32 (vjw9): switched ps1 in
        # _ps1_preserved / is_emission_canonical / same_emission_orbit /
        # count_dedup_prefilter to the right action — all three dedup
        # strategies now agree at 18 (the qgraf-validated count).
        m = qed_model()
        @test count_dedup_burnside(m,  [:e,:e], [:mu,:mu]; loops=1) == 18
        @test count_dedup_canonical(m, [:e,:e], [:mu,:mu]; loops=1) == 18
        @test count_dedup_prefilter(m, [:e,:e], [:mu,:mu]; loops=1) == 18
    end

    @testset "AUDITION VERDICT: all three dedup strategies agree on QED ee→ee" begin
        # QED1 ee→ee tree: legacy = 2 (s + t channels).
        # Session 32 (bead vjw9) closed the over-dedup bug in (B) and (C):
        # the ps1 action was incorrectly LEFT-action `g[ps1[i]]`, which
        # relabelled physical legs instead of slots. With the fix all
        # three strategies produce 2.
        m = qed1_model()
        @test count_dedup_burnside(m,  [:e,:e], [:e,:e]; loops=0) == 2
        @test count_dedup_canonical(m, [:e,:e], [:e,:e]; loops=0) == 2
        @test count_dedup_prefilter(m, [:e,:e], [:e,:e]; loops=0) == 2
    end

end
