#  Phase 17b: full battery — count_diagrams_qg21 vs legacy on every test
#  in test_diagram_gen.jl.  Same cases that appear in the in-tree golden
#  master tests; if Strategy C is ready to ship, both must agree.

using Test
using Feynfeld
using Feynfeld.QgrafPort: count_diagrams_qg21

@testset "Phase 17b: qg21 path vs legacy on diagram_gen battery" begin

    # ─── φ³ tree-level ───────────────────────────────────────────────
    @testset "φ³ tree-level (5 cases)" begin
        m = phi3_model()
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi]; loops=0) == 1
        @test count_diagrams_qg21(m, [:phi, :phi], [:phi, :phi]; loops=0) == 3
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi, :phi]; loops=0) == 3
        @test count_diagrams_qg21(m, [:phi, :phi], [:phi, :phi, :phi]; loops=0) == 15
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi, :phi, :phi]; loops=0) == 15
    end

    @testset "φ³ 1-loop (4 cases)" begin
        m = phi3_model()
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi]; loops=1) == 7
        @test count_diagrams_qg21(m, [:phi, :phi], [:phi, :phi]; loops=1) == 39
        @test count_diagrams_qg21(m, [:phi, :phi], [:phi, :phi]; loops=1, onepi=true) == 3
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi, :phi]; loops=1) == 39
    end

    @testset "φ³ 2-loop — THE 465 regression case (fixed Phase 12e)" begin
        # BUG 2 (Session 24): Julia qg21 emitted 52 canonical topologies for
        # the 6-deg-3 partition; qgraf emits 50.  Root cause: Julia's
        # step_c_enumerate! lacked qgraf's post-fill permutation canonicality
        # check (qg21:13156-13291).  Step C's cross-row/col checks
        # (qg21:12911-12946) are necessary but not sufficient.  Fix:
        # is_canonical_qgraf! in canonical.jl + call from emit path.
        m = phi3_model()
        @test count_diagrams_qg21(m, [:phi, :phi], [:phi, :phi]; loops=2) == 465
    end

    @testset "φ³ 2-loop — φ→φφ (single-incoming case agrees)" begin
        m = phi3_model()
        @test count_diagrams_qg21(m, [:phi], [:phi, :phi]; loops=2) == 58
    end

    @testset "QED 1-gen tree-level (5 cases)" begin
        m = qed1_model()
        @test count_diagrams_qg21(m, [:e], [:e, :gamma]; loops=0) == 1
        @test count_diagrams_qg21(m, [:gamma], [:e, :e]; loops=0) == 1
        @test count_diagrams_qg21(m, [:e, :e], [:e, :e]; loops=0) == 2     # s + t
        @test count_diagrams_qg21(m, [:e, :gamma], [:e, :gamma]; loops=0) == 2
        @test count_diagrams_qg21(m, [:e, :e], [:gamma, :gamma]; loops=0) == 2
    end

    @testset "QED 2-gen (e+μ) tree (4 cases)" begin
        m = qed_model()
        @test count_diagrams_qg21(m, [:e, :e], [:mu, :mu]; loops=0) == 1
        @test count_diagrams_qg21(m, [:e, :mu], [:e, :mu]; loops=0) == 1
        @test count_diagrams_qg21(m, [:mu, :mu], [:e, :e]; loops=0) == 1
        @test count_diagrams_qg21(m, [:e, :mu], [:mu, :e]; loops=0) == 1
    end

    @testset "QCD tree-level" begin
        m = qcd_model()
        @test count_diagrams_qg21(m, [:q, :q], [:g, :g]; loops=0) == 3
        @test count_diagrams_qg21(m, [:g, :g], [:g, :g]; loops=0) == 4   # with gggg
        @test count_diagrams_qg21(m, [:q, :g], [:q, :g]; loops=0) == 3
        @test count_diagrams_qg21(m, [:g, :g], [:q, :q]; loops=0) == 3
    end

end
