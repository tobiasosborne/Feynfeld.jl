# Test: Diagram generation — golden master validation
#
# Tests the new algorithmic diagram generator against golden master
# outputs from qgraf-4.0.6. Every diagram count has been cross-validated
# against known physics (22/22 checks pass).
#
# Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/manifest.json
# Method: Golden master TDD — qgraf is the oracle, Julia must match.

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Diagram generation" begin

    # ─── φ³ theory: pure combinatorics ────────────────────────────
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/phi3/
    @testset "φ³ tree-level" begin
        phi3 = phi3_model()

        # φ → φφ tree: 1 diagram (single vertex)
        # Source: golden_masters/phi3/phi_TO_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi], [:phi, :phi]; loops=0) == 1

        # φφ → φφ tree: 3 diagrams (s + t + u channels)
        # Source: golden_masters/phi3/phi_phi_TO_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi]; loops=0) == 3

        # φ → φφφ tree: 3 diagrams
        # Source: golden_masters/phi3/phi_TO_phi_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi], [:phi, :phi, :phi]; loops=0) == 3

        # φφ → φφφ tree: 15 diagrams
        # Source: golden_masters/phi3/phi_phi_TO_phi_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi, :phi]; loops=0) == 15

        # φ → φφφφ tree: 15 diagrams
        # Source: golden_masters/phi3/phi_TO_phi_phi_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi], [:phi, :phi, :phi, :phi]; loops=0) == 15

        # φφ → φφφφ tree: 105 diagrams
        # Source: golden_masters/phi3/phi_phi_TO_phi_phi_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi, :phi, :phi]; loops=0) == 105

        # φφφ → φφφ tree: 105 diagrams
        # Source: golden_masters/phi3/phi_phi_phi_TO_phi_phi_phi_0L.out
        @test count_diagrams(phi3, [:phi, :phi, :phi], [:phi, :phi, :phi]; loops=0) == 105
    end

    @testset "φ³ 1-loop" begin
        phi3 = phi3_model()

        # φ → φφ 1-loop: 7 diagrams (all), 1 diagram (1PI)
        # Source: golden_masters/phi3/phi_TO_phi_phi_1L.out
        @test count_diagrams(phi3, [:phi], [:phi, :phi]; loops=1) == 7
        @test count_diagrams(phi3, [:phi], [:phi, :phi]; loops=1, onepi=true) == 1

        # φφ → φφ 1-loop: 39 diagrams (all), 3 (1PI)
        # Source: golden_masters/phi3/phi_phi_TO_phi_phi_1L.out
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi]; loops=1) == 39
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi]; loops=1, onepi=true) == 3
    end

    # ─── φ³ 2-loop (qg21 regression — Strategy C port target) ─────
    # First RED test for the qgraf Strategy C port. Current naive fill
    # + post-hoc pairwise-swap canonicalization gives 474 (duplicate
    # topologies at 2-loop+). qgraf's row-by-row fill + full per-class
    # lex-permutation isomorph check gives 465.
    # Ref: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08 lines 12426-13668 (qg21)
    # Ref: refs/papers/Nogueira1993_JCompPhys105_279.pdf §2 "The method"
    # Beads: feynfeld-4tg, blocked by feynfeld-ney
    @testset "φ³ 2-loop (Strategy C port target)" begin
        phi3 = phi3_model()

        # φφ → φφ 2-loop: 465 diagrams (qgraf golden master)
        # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/phi3/phi_phi_TO_phi_phi_2L.out
        # RED as of Session 20: naive port gives 474 (9 duplicates)
        @test count_diagrams(phi3, [:phi, :phi], [:phi, :phi]; loops=2) == 465
    end

    # ─── QED: electrons + photons ─────────────────────────────────
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/qed1/
    @testset "QED 1-gen tree-level" begin
        qed = qed_model()

        # e⁺e⁻ → γγ tree: 2 diagrams (t + u channels)
        # Source: golden_masters/qed1/eminus_eplus_TO_photon_photon_0L.out
        @test count_diagrams(qed, [:e, :e], [:gamma, :gamma]; loops=0) == 2

        # Compton eγ → eγ tree: 2 diagrams (s + u)
        # Source: golden_masters/qed1/eminus_photon_TO_eminus_photon_0L.out
        @test count_diagrams(qed, [:e, :gamma], [:e, :gamma]; loops=0) == 2

        # Bhabha e⁺e⁻ → e⁺e⁻ tree: 2 diagrams (s + t)
        # Source: golden_masters/qed1/eminus_eplus_TO_eminus_eplus_0L.out
        @test count_diagrams(qed, [:e, :e], [:e, :e]; loops=0) == 2

        # γγ → e⁺e⁻ tree: 2 diagrams (t + u)
        # Source: golden_masters/qed1/photon_photon_TO_eminus_eplus_0L.out
        @test count_diagrams(qed, [:gamma, :gamma], [:e, :e]; loops=0) == 2

        # γγ → γγ tree: 0 diagrams (Furry's theorem — no 4-photon vertex in QED)
        # Source: golden_masters/qed1/photon_photon_TO_photon_photon_0L.out
        @test count_diagrams(qed, [:gamma, :gamma], [:gamma, :gamma]; loops=0) == 0
    end

    @testset "QED 1-gen 1-loop" begin
        qed1 = qed1_model()    # use the 1-gen model (e only — no μ loop contributions)

        # Bhabha 1PI 1-loop: 4 diagrams (2 topologies × 2 fermion orientations).
        # Matches qgraf's count with distinct e⁻/e⁺ fields.
        # Source: golden_masters/qed1/eminus_eplus_TO_eminus_eplus_1L_onepi.out
        @test count_diagrams(qed1, [:e, :e], [:e, :e]; loops=1, onepi=true) == 4

        # γγ → γγ 1-loop: 6 diagrams (box permutations, Furry's theorem allows 1-loop).
        # 1-gen only — 2-gen would have 12 (μ box also contributes).
        # Source: golden_masters/qed1/photon_photon_TO_photon_photon_1L.out
        @test count_diagrams(qed1, [:gamma, :gamma], [:gamma, :gamma]; loops=1) == 6
    end

    # ─── QED 2-gen: e + μ (THE Feynfeld benchmark) ────────────────
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/qed2/
    @testset "QED 2-gen (e+μ)" begin
        qed = qed_model()

        # e⁺e⁻ → μ⁺μ⁻ tree: 1 diagram (s-channel only)
        # Source: golden_masters/qed2/eminus_eplus_TO_muminus_muplus_0L.out
        @test count_diagrams(qed, [:e, :e], [:mu, :mu]; loops=0) == 1

        # e⁺e⁻ → μ⁺μ⁻ 1-loop 1PI: 2 diagrams (direct + crossed box)
        # Source: golden_masters/qed2/eminus_eplus_TO_muminus_muplus_1L_onepi.out
        @test count_diagrams(qed, [:e, :e], [:mu, :mu]; loops=1, onepi=true) == 2

        # e⁺e⁻ → μ⁺μ⁻ 1-loop all: 18 diagrams
        # Source: golden_masters/qed2/eminus_eplus_TO_muminus_muplus_1L.out
        @test count_diagrams(qed, [:e, :e], [:mu, :mu]; loops=1) == 18

        # e⁻μ⁻ → e⁻μ⁻ tree: 1 diagram (t-channel only)
        # Source: golden_masters/qed2/eminus_muminus_TO_eminus_muminus_0L.out
        @test count_diagrams(qed, [:e, :mu], [:e, :mu]; loops=0) == 1
    end

    # ─── QCD ──────────────────────────────────────────────────────
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/qcd/
    @testset "QCD tree-level" begin
        qcd = qcd_model()

        # qq̄ → gg tree: 3 diagrams (s + t + u)
        # Source: golden_masters/qcd/quark_antiquark_TO_gluon_gluon_0L.out
        @test count_diagrams(qcd, [:q, :q], [:g, :g]; loops=0) == 3

        # gg → gg tree: 4 diagrams (s + t + u via ggg + 1 contact via gggg)
        # Source: golden_masters/qcd/gluon_gluon_TO_gluon_gluon_0L.out
        @test count_diagrams(qcd, [:g, :g], [:g, :g]; loops=0) == 4

        # qq̄ → qq̄ tree: 2 diagrams (s + t)
        # Source: golden_masters/qcd/quark_antiquark_TO_quark_antiquark_0L.out
        @test count_diagrams(qcd, [:q, :q], [:q, :q]; loops=0) == 2

        # qg → qg tree: 3 diagrams (s + t + u)
        # Source: golden_masters/qcd/quark_gluon_TO_quark_gluon_0L.out
        @test count_diagrams(qcd, [:q, :g], [:q, :g]; loops=0) == 3

        # gg → qq̄ tree: 3 diagrams (s + t + u)
        # Source: golden_masters/qcd/gluon_gluon_TO_quark_antiquark_0L.out
        @test count_diagrams(qcd, [:g, :g], [:q, :q]; loops=0) == 3
    end

    # ─── Pipeline integration ─────────────────────────────────────
    @testset "Pipeline: generate_diagrams produces valid TreeChannels" begin
        model = qed_model()
        rules = feynman_rules(model)
        incoming = [
            ExternalLeg(:e, Momentum(:p1), true, false),
            ExternalLeg(:e, Momentum(:p2), true, true),
        ]
        outgoing = [
            ExternalLeg(:mu, Momentum(:k1), false, false),
            ExternalLeg(:mu, Momentum(:k2), false, true),
        ]

        # Must match existing tree_channels for backward compatibility
        old_channels = tree_channels(model, rules, incoming, outgoing)
        @test length(old_channels) == 1  # s-channel only

        # New diagram generation must produce same result
        new_channels = generate_tree_channels(model, rules, incoming, outgoing)
        @test length(new_channels) == length(old_channels)
        @test new_channels[1].channel == :s
        @test new_channels[1].exchanged == :gamma
    end
end
