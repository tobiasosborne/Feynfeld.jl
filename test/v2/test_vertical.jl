# Vertical tracer bullet: e+e- → μ+μ- through ALL layers
#
# Layer 1 (Model)    → define QED
# Layer 2 (Rules)    → extract Feynman rules
# Layer 3 (Diagrams) → generate tree diagram
# Layer 4 (Algebra)  → spin sum, trace, contract → |M|²
# Layer 6 (Evaluate) → cross-section σ = 4πα²/(3s)
#
# Validates: P&S Eq (5.10): |M|² = 8(t²+u²)
#            P&S Eq (5.12): σ = 4πα²/(3s)

using Test
using Feynfeld

@testset "Vertical: e+e- → μ+μ-" begin

    # ======== Layer 1: Model ========
    @testset "Layer 1: QED Model" begin
        model = qed_model()
        @test model isa QEDModel
        @test model isa AbstractModel
        @test model_name(model) == :QED
        @test length(model_fields(model)) == 3

        # Field properties via dispatch
        e = get_field(model, :e)
        @test e isa Field{Fermion}
        @test species(e) isa Fermion
        @test mass_trait(e) == Massive()
        @test charge_trait(e) == Charged()

        γ = get_field(model, :gamma)
        @test γ isa Field{Boson}
        @test mass_trait(γ) == Massless()

        # Gauge group
        @test gauge_groups(model) == [U1()]

        # Derived interface
        @test length(fermion_fields(model)) == 2
        @test length(boson_fields(model)) == 1
    end

    # ======== Layer 2: Rules ========
    @testset "Layer 2: Feynman Rules" begin
        model = qed_model()
        rules = feynman_rules(model)
        @test rules isa FeynmanRules

        # Should have vertices for e-e-γ and μ-μ-γ
        @test length(rules.vertices) == 2

        # Callable: look up vertex
        v = rules((:e, :e, :gamma))
        @test v isa VertexRule
        @test v.coupling == :e

        # Vertex factor produces DiracExpr
        mu = LorentzIndex(:mu)
        vf = vertex_factor(rules, (:e, :e, :gamma), mu)
        @test vf isa DiracExpr

        # Propagator dispatch on species
        p = Momentum(:p)
        prop_e = propagator_num(Fermion(), p, 0//1)
        @test prop_e isa DiracExpr

        prop_γ = propagator_num(Boson(), mu, LorentzIndex(:nu))
        @test prop_γ isa AlgSum
    end

    # ======== Layer 3: Channels ========
    @testset "Layer 3: Tree Channels" begin
        model = qed_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [
            ExternalLeg(:e, p1, true, false),   # e-
            ExternalLeg(:e, p2, true, true),     # e+ (antiparticle)
        ]
        outgoing = [
            ExternalLeg(:mu, k1, false, false),  # μ-
            ExternalLeg(:mu, k2, false, true),    # μ+
        ]

        channels = tree_channels(model, rules, incoming, outgoing)
        @test length(channels) == 1
        @test channels[1].channel == :s
        @test channels[1].exchanged == :gamma

        # Build amplitude
        chain_L, chain_R = build_amplitude(channels[1], rules, model)
        @test chain_L isa DiracExpr
        @test chain_R isa DiracExpr
        @test length(chain_L.terms) == 1  # single vertex term (QED)
        @test length(chain_L.terms[1][2].elements) == 3  # vbar, gamma, u
    end

    # ======== Layer 4: Algebra (validated before, run again) ========
    @testset "Layer 4: Spin sum → trace → contract" begin
        model = qed_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
        outgoing = [ExternalLeg(:mu, k1, false, false), ExternalLeg(:mu, k2, false, true)]

        channels = tree_channels(model, rules, incoming, outgoing)
        chain_L, chain_R = build_amplitude(channels[1], rules, model)

        # Spin-averaged |M|²
        m_sq = spin_sum_amplitude_squared(chain_L, chain_R)
        @test m_sq isa AlgSum
        @test !iszero(m_sq)

        # Contract Lorentz indices
        contracted = contract(m_sq)
        @test contracted isa AlgSum

        # Expand scalar products
        expanded = expand_scalar_product(contracted)
        @test expanded isa AlgSum
    end

    # ======== Layer 6: Cross-section ========
    @testset "Layer 6: Evaluate → P&S (5.10) and (5.12)" begin
        # Full pipeline via solve_tree
        model = qed_model()
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
        outgoing = [ExternalLeg(:mu, k1, false, false), ExternalLeg(:mu, k2, false, true)]

        prob = CrossSectionProblem(model, incoming, outgoing, 100)
        result = solve_tree(prob)

        # Evaluate at specific kinematics: s=10, t=-3, u=-7
        man = Mandelstam(10, -3//1, -7//1)
        m_sq = evaluate_m_squared(result.amplitude_squared, man)

        # P&S (5.10): |M|² = 8(t² + u²) (coupling stripped)
        expected = 8 * ((-3)^2 + (-7)^2)
        @test m_sq == expected  # 8(9+49) = 464

        # P&S (5.12): σ_total = 4πα²/(3s) for massless e+e- → μ+μ-
        α = 1/137.036
        s = 91.2^2  # GeV² (Z-pole energy)
        σ = sigma_total_tree_ee_mumu(s; alpha=α)
        # Convert to nanobarns: 1 GeV⁻² = 0.3894e6 nb
        σ_nb = σ * 0.3894e6
        # QED-only (no Z): σ ≈ 0.0104 nb at √s = 91.2 GeV
        # (The experimental ~87 nb is from Z resonance, not pure QED)
        @test 0.005 < σ_nb < 0.02

        # At √s = 10 GeV: σ ≈ 0.87 nb = 869 pb
        σ_10 = sigma_total_tree_ee_mumu(100.0; alpha=α)
        σ_10_nb = σ_10 * 0.3894e6
        @test 0.5 < σ_10_nb < 1.5

        println("  Vertical tracer bullet results:")
        println("  |M|² at s=10,t=-3,u=-7: $m_sq (expected $expected)")
        println("  σ(√s=91.2 GeV) = $(round(σ_nb, digits=3)) nb")
        println("  σ(√s=10 GeV)   = $(round(σ_10_nb*1000, digits=1)) pb")
    end
end
