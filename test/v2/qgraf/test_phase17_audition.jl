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
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: count_dedup_burnside, count_dedup_canonical,
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

    @testset "QED ee→μμ 1L (2-gen) — Burnside fixed (Phase 12d)" begin
        # Phase 12d (qgen.jl): _qgen_recurse now enumerates all positional
        # perms of remaining slots + applies self-loop & multi-edge filters
        # (qgen:13921-13954).  Was A=16, B=17, C=17; now A=18, B=19, C=19.
        # (A) Burnside is now correct.  (B) and (C) over-count because
        # their canonicality compare-under-autos sees the new orbit reps as
        # all-distinct rather than orbit-equivalent — the known canonical-rep
        # bug from the audition VERDICT (HANDOFF.md), unaffected by Phase 12d.
        m = qed_model()
        @test         count_dedup_burnside(m,  [:e,:e], [:mu,:mu]; loops=1) == 18
        @test_broken  count_dedup_canonical(m, [:e,:e], [:mu,:mu]; loops=1) == 18
        @test_broken  count_dedup_prefilter(m, [:e,:e], [:mu,:mu]; loops=1) == 18
    end

    @testset "AUDITION VERDICT: (A) Burnside agrees with legacy on QED ee→ee" begin
        # QED1 ee→ee tree: legacy = 2 (s + t channels).
        # (A) Burnside    = 2  ✓
        # (B) canonical   = 1  ✗ over-dedups (s ≡ t under in↔out auto, wrong)
        # (C) pre-filter  = 1  ✗ same bug
        m = qed1_model()
        @test count_dedup_burnside(m, [:e, :e], [:e, :e]; loops=0) == 2
        @test_broken count_dedup_canonical(m, [:e, :e], [:e, :e]; loops=0) == 2
        @test_broken count_dedup_prefilter(m, [:e, :e], [:e, :e]; loops=0) == 2
    end

end
