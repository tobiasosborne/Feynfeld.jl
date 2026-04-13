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

    @testset "[broken] QED ee→μμ 1L 2gen — qgen flavor-loop bug (Phase 12d)" begin
        # Legacy = 18; qgen returns ≤17.  See audition commit for diagnosis.
        @test_broken count_diagrams_qg21(qed_model(), [:e, :e], [:mu, :mu]; loops=1) ==
                     count_diagrams(qed_model(),     [:e, :e], [:mu, :mu]; loops=1)
    end

end
