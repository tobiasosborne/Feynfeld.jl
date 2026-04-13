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

    @testset "QED ee→μμ 1L (2-gen) — discrepancy across all 3" begin
        # Audition reveals a real bug: legacy 18 vs (A=16, B=17, C=17).
        # Diagnosis: my qgen multiset-matching seems to miss flavor-loop
        # variants (e-loop AND μ-loop) for some perms.  Investigate before
        # the pipeline swap (Phase 17b) — this is the kind of bug the
        # audition was designed to catch.
        m = qed_model()
        @test_broken count_dedup_burnside(m,  [:e,:e], [:mu,:mu]; loops=1) == 18
        @test_broken count_dedup_canonical(m, [:e,:e], [:mu,:mu]; loops=1) == 18
        @test_broken count_dedup_prefilter(m, [:e,:e], [:mu,:mu]; loops=1) == 18
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
