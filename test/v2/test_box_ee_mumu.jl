# Spiral 10: 1-loop box diagrams for e+e- → μ+μ- (QED, massless)
#
# Ground truth: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/ElAel-MuAmu.m
# Cross-validated against arXiv:hep-ph/0010075, Eq. 2.32 (Ellis, Giele, Zanderighi)
#
# Direct box: D₀(0,0,0,0,s,t, 0,0,0,0)
# Crossed box: D₀(0,0,0,0,s,u, 0,0,0,0)
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 4.4

using Test
using Feynfeld

# ---- Setup: ee→μμ process ----
const p1 = Momentum(:p1)
const p2 = Momentum(:p2)
const k1 = Momentum(:k1)
const k2 = Momentum(:k2)

const incoming_legs = [
    ExternalLeg(:e, p1, true, false),   # e⁻ incoming
    ExternalLeg(:e, p2, true, true),    # e⁺ incoming (antiparticle)
]
const outgoing_legs = [
    ExternalLeg(:mu, k1, false, false),  # μ⁻ outgoing
    ExternalLeg(:mu, k2, false, true),   # μ⁺ outgoing (antiparticle)
]

@testset "Spiral 10: 1-loop box infrastructure" begin

    # ==== Phase A: LoopChannel types and enumeration ====
    @testset "Phase A: box_channels enumeration" begin
        model = qed_model()
        rules = feynman_rules(model)
        channels = box_channels(model, rules, incoming_legs, outgoing_legs)

        @test length(channels) == 2
        topos = Set(ch.topology for ch in channels)
        @test :direct_box ∈ topos
        @test :crossed_box ∈ topos

        # Both channels should have loop momentum :q
        @test all(ch.loop_momentum == :q for ch in channels)

        # Internal fields: 2 photon + 2 fermion
        for ch in channels
            gamma_count = count(==(:gamma), ch.internal_fields)
            @test gamma_count == 2
        end
    end

    # ==== Phase B: Loop amplitude construction ====
    @testset "Phase B: build_loop_box_amplitude" begin
        model = qed_model()
        rules = feynman_rules(model)
        channels = box_channels(model, rules, incoming_legs, outgoing_legs)

        for ch in channels
            chain_e, chain_mu, denoms = build_loop_box_amplitude(ch, rules, model)

            # Both chains should be single-term DiracExpr
            @test length(chain_e.terms) == 1
            @test length(chain_mu.terms) == 1

            # Each chain should have 5 elements: bar_sp, γ, GS(q+p), γ, plain_sp
            _, dc_e = chain_e.terms[1]
            _, dc_mu = chain_mu.terms[1]
            @test length(dc_e.elements) == 5
            @test length(dc_mu.elements) == 5

            # First and last elements should be spinors
            @test dc_e.elements[1] isa Spinor
            @test dc_e.elements[5] isa Spinor
            @test dc_mu.elements[1] isa Spinor
            @test dc_mu.elements[5] isa Spinor

            # Middle element should contain loop momentum (MomSumSlot)
            @test dc_e.elements[3] isa DiracGamma{MomSumSlot}
            @test dc_mu.elements[3] isa DiracGamma{MomSumSlot}

            # BoxDenominators: all masses zero (massless QED)
            @test denoms.masses == (0//1, 0//1, 0//1, 0//1)
        end

        # Direct box: p₁₃ = t, crossed: p₁₃ = u
        ch_direct = first(ch for ch in channels if ch.topology == :direct_box)
        ch_crossed = first(ch for ch in channels if ch.topology == :crossed_box)
        _, _, d_direct = build_loop_box_amplitude(ch_direct, rules, model)
        _, _, d_crossed = build_loop_box_amplitude(ch_crossed, rules, model)
        @test d_direct.kinematic_inv == :t
        @test d_crossed.kinematic_inv == :u
    end

    # ==== Phase C: Tree × loop interference spin sum ====
    @testset "Phase C: spin_sum_tree_loop_interference" begin
        model = qed_model()
        rules = feynman_rules(model)

        # Build tree amplitude (s-channel photon exchange)
        tree_channels_list = tree_channels(model, rules, incoming_legs, outgoing_legs)
        @test !isempty(tree_channels_list)
        tree_ch = first(tree_channels_list)
        tree_amp_L, tree_amp_R = build_amplitude(tree_ch, rules, model)

        # Build loop amplitude (direct box)
        loop_channels_list = box_channels(model, rules, incoming_legs, outgoing_legs)
        ch_direct = first(ch for ch in loop_channels_list if ch.topology == :direct_box)
        loop_e, loop_mu, denoms = build_loop_box_amplitude(ch_direct, rules, model)

        # Compute interference: Σ_spins M_tree* × M_box
        interference = spin_sum_tree_loop_interference(
            (tree_amp_L, tree_amp_R),
            (loop_e, loop_mu))

        # The result should be a non-empty AlgSum
        @test interference isa AlgSum
        @test !isempty(interference.terms)

        # Contract Lorentz indices
        contracted = contract(interference)
        @test contracted isa AlgSum

        # Expand scalar products (including loop momentum q)
        expanded = expand_scalar_product(contracted)
        @test expanded isa AlgSum

        # After expansion, all Lorentz indices should be contracted.
        # Only ScalarProduct (Pair{Momentum,Momentum}) should remain.
        for (fk, _) in expanded.terms
            for f in fk.factors
                @test f isa ScalarProduct
            end
        end

        # The result should contain SP(q, ...) terms (loop momentum dependent)
        function _has_loop_momentum(terms)
            for (fk, _) in terms
                for f in fk.factors
                    if f isa ScalarProduct
                        (f.a.name == :q || f.b.name == :q) && return true
                    end
                end
            end
            false
        end
        @test _has_loop_momentum(expanded.terms)

        # Count the number of terms (rough sanity check)
        # Trace of 6 gammas × 6 gammas produces O(100) terms
        @test 10 < length(expanded.terms) < 500
    end
end
