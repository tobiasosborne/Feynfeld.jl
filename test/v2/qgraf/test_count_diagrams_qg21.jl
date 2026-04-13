#  Phase 17b: count_diagrams_qg21 — new Strategy C entry point.
#
#  Wraps qg21 + qg10 + qgen + (A) Burnside dedup.  Same signature as
#  legacy count_diagrams; should agree on every case where neither side
#  has a known bug.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: count_diagrams_qg21

@testset "Phase 17b: count_diagrams_qg21 vs legacy" begin

    @testset "phi3 tree φ→φφ" begin
        @test count_diagrams_qg21(phi3_model(), [:phi], [:phi, :phi]; loops=0) ==
              count_diagrams(phi3_model(),     [:phi], [:phi, :phi]; loops=0)
    end

    @testset "phi3 tree φφ→φφ" begin
        @test count_diagrams_qg21(phi3_model(), [:phi, :phi], [:phi, :phi]; loops=0) ==
              count_diagrams(phi3_model(),     [:phi, :phi], [:phi, :phi]; loops=0)
    end

    @testset "phi3 tree φ→φφφ" begin
        @test count_diagrams_qg21(phi3_model(), [:phi], [:phi, :phi, :phi]; loops=0) ==
              count_diagrams(phi3_model(),     [:phi], [:phi, :phi, :phi]; loops=0)
    end

    @testset "phi3 1L φ→φφ" begin
        @test count_diagrams_qg21(phi3_model(), [:phi], [:phi, :phi]; loops=1) ==
              count_diagrams(phi3_model(),     [:phi], [:phi, :phi]; loops=1)
    end

    @testset "phi3 1L φφ→φφ" begin
        @test count_diagrams_qg21(phi3_model(), [:phi, :phi], [:phi, :phi]; loops=1) ==
              count_diagrams(phi3_model(),     [:phi, :phi], [:phi, :phi]; loops=1)
    end

    @testset "QED ee→μμ tree" begin
        @test count_diagrams_qg21(qed_model(), [:e, :e], [:mu, :mu]; loops=0) ==
              count_diagrams(qed_model(),     [:e, :e], [:mu, :mu]; loops=0)
    end

    @testset "QED1 ee→ee tree (s + t)" begin
        @test count_diagrams_qg21(qed1_model(), [:e, :e], [:e, :e]; loops=0) ==
              count_diagrams(qed1_model(),     [:e, :e], [:e, :e]; loops=0)
    end

    @testset "QED1 eγ→eγ tree" begin
        @test count_diagrams_qg21(qed1_model(), [:e, :gamma], [:e, :gamma]; loops=0) ==
              count_diagrams(qed1_model(),     [:e, :gamma], [:e, :gamma]; loops=0)
    end

    @testset "onepi filter: phi3 1L φφ→φφ onepi → 3" begin
        # Legacy count_diagrams supports onepi flag too.
        @test count_diagrams_qg21(phi3_model(), [:phi, :phi], [:phi, :phi];
                                    loops=1, onepi=true) ==
              count_diagrams(phi3_model(),     [:phi, :phi], [:phi, :phi];
                                    loops=1, onepi=true)
    end

    @testset "QED ee→μμ 1L 2gen — flavor-loop bug fix (Phase 12d)" begin
        # Was @test_broken — Burnside returned 16 vs legacy 18.
        # Fixed by enumerating all positional perms of remaining slots in
        # _qgen_recurse + applying self-loop & multi-edge filters
        # (qgen:13921-13954).  See grind/ for the diagnostic infrastructure.
        @test count_diagrams_qg21(qed_model(), [:e, :e], [:mu, :mu]; loops=1) ==
              count_diagrams(qed_model(),     [:e, :e], [:mu, :mu]; loops=1)
    end

end
