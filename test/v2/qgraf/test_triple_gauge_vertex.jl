# Triple gauge-boson vertex in the qgraf amplitude pipeline (bead
# feynfeld-asp9). Wires the ggg vertex into build_vertices and assembles
# the boson sub-amplitude (off-fermion-line vertex factors) in
# emission_to_amplitude, so qq̄→gg flows through solve_tree_pipeline.
#
# Ground truth:
#   - Triple-gauge Lorentz structure: refs/papers/PeskinSchroeder1995.djvu
#     Eq. (16.10), implemented by `triple_gauge_vertex` (gauge_vertex.jl).
#   - qq̄→gg channel decomposition: test/v2/test_qqbar_gg.jl, itself
#     validated against refs/FeynCalc/.../QCD/Tree/.../QQbar-GlGl.m.
#
# Three tree diagrams: s-channel (triple-gluon vertex), t-channel and
# u-channel (quark exchange). The t/u channels are pure multi-vertex
# fermion lines and already flowed through the pipeline (Phase 18b-3 +
# 18b-4); this bead adds the s-channel's off-fermion-line ggg vertex.
#
# Per-channel diagonal traces are checked symbolically against handbuilt
# references built from the cited `triple_gauge_vertex` + `dirac_trace`,
# with Feynman-gauge polarisation sums (matching combine_m_squared_burnside).
# Full qq̄→gg ≡ FeynCalc additionally needs colour algebra in the pipeline
# (bead feynfeld-yewo) and is the remit of feynfeld-4xrh.

using Test

using Feynfeld
using Feynfeld.QgrafPort: AmplitudeBundle, emission_to_amplitude,
                            enumerate_topology_automorphisms,
                            is_emission_canonical, _foreach_emission,
                            _pair_trace

@testset "Triple gauge vertex: qq̄→gg through the pipeline" begin

    p1 = Momentum(:p1); p2 = Momentum(:p2)
    k1 = Momentum(:k1); k2 = Momentum(:k2)

    # Collect the canonical qq̄→gg emissions — mirrors solve_tree_pipeline's
    # own emission loop, but keeps the per-bundle structure for inspection.
    function qqgg_bundles()
        model = qcd_model(m_q=:zero)
        legs = [ExternalLeg(:q, p1, true, false), ExternalLeg(:q, p2, true, true),
                ExternalLeg(:g, k1, false, false), ExternalLeg(:g, k2, false, false)]
        ext_exp = Feynfeld._qgraf_ext_labels(model, legs)
        pm = [l.momentum for l in legs]
        pa = Bool[l.antiparticle for l in legs]
        bundles = AmplitudeBundle[]
        _foreach_emission(model, [:q, :q_bar], [:g, :g];
                          loops=0, ext_exp=ext_exp) do st, lb, ps1, pmap
            is_emission_canonical(st, lb, enumerate_topology_automorphisms(st),
                                  ps1, pmap) || return
            push!(bundles, emission_to_amplitude(st, lb, ps1, pmap, model;
                    physical_moms=pm, n_inco=2, phys_anti=pa))
        end
        bundles
    end
    bundles = qqgg_bundles()

    @testset "qq̄→gg emits 3 canonical orbits (s, t, u)" begin
        @test length(bundles) == 3
    end

    s_b, t_b, u_b = bundles   # emission order: s-channel, t-channel, u-channel

    eps3  = LorentzIndex(:eps_3,  DimD()); eps3c = LorentzIndex(:eps_3_,  DimD())
    eps4  = LorentzIndex(:eps_4,  DimD()); eps4c = LorentzIndex(:eps_4_,  DimD())
    mu5   = LorentzIndex(:mu_l_5, DimD()); mu5c  = LorentzIndex(:mu_l_5_, DimD())
    γ(idx) = DiracGamma(LISlot(idx))

    @testset "s-channel: ggg vertex assembled into boson_factor" begin
        # s-channel = qqg vertex on the fermion line (γ^{:mu_l_5}) plus the
        # ggg vertex OFF the fermion line, captured in boson_factor.
        @test length(s_b.line_chains) == 1
        @test s_b.boson_pols == [eps3, eps4]
        # boson_factor = the P&S 16.10 triple-gauge vertex, external-gluon
        # indices canonicalised to :eps_3/:eps_4, internal index :mu_l_5.
        # Momenta all-outgoing from v6: external gluons k1,k2 leave the
        # diagram → +k1,+k2; the internal edge carries p1+p2 INTO v6 → -(p1+p2).
        negq = MomentumSum([(-1//1, p1), (-1//1, p2)])
        @test s_b.boson_factor == triple_gauge_vertex(eps3, eps4, mu5, k1, k2, negq)
        # internal-gluon propagator denominator (p1+p2)².
        q = MomentumSum([(1//1, p1), (1//1, p2)])
        @test s_b.denoms == [alg(pair(q, q))]
    end

    @testset "t/u-channels: pure fermion line, trivial boson_factor" begin
        @test t_b.boson_factor == alg(1)
        @test u_b.boson_factor == alg(1)
        @test isempty(s_b.line_chains) == false   # all three have a quark line
        @test length(t_b.line_chains) == 1 && length(u_b.line_chains) == 1
    end

    @testset "per-bundle coupling = g_s² (Phase 18b-7)" begin
        # Each diagram has 2 g_s vertices: s = qqg + ggg, t/u = 2× qqg.
        # `VertexRule.coupling_power = 1` for each, so total = 2.
        for b in (s_b, t_b, u_b)
            @test b.coupling == coupling_alg(:g_s, 2)
        end
    end

    # Per-channel diagonal |M_X|² ≡ handbuilt reference. Handbuilt traces
    # are built from the cited triple_gauge_vertex (s) / explicit quark-
    # exchange γ-chains (t,u) plus dirac_trace, with Feynman-gauge ε–ε*
    # sums — exactly the structure combine_m_squared_burnside assembles.
    # Phase 18b-7 (feynfeld-5d1k): every qq̄→gg diagram has 2 g_s vertices
    # (s: qqg + ggg; t/u: 2× qqg), so each bundle carries `g_s²` and `_pair_trace`
    # multiplies in `g_s²·g_s² = g_s⁴` for every (i,j) pair.
    polsum()  = polarization_sum(eps3, eps3c) * polarization_sum(eps4, eps4c)
    finish(x) = expand_scalar_product(contract(x * polsum()))
    gs4 = coupling_alg(:g_s, 4)

    @testset "s-channel diagonal |M_s|² ≡ handbuilt D_ss" begin
        negq = MomentumSum([(-1//1, p1), (-1//1, p2)])
        # quark line ū(p2) γ^{:mu_l_5} u(p1): Tr[p̸1 γ^{:mu_l_5_} p̸2 γ^{:mu_l_5}].
        tr = dirac_trace(DiracGamma[GS(p1), γ(mu5c), GS(p2), γ(mu5)])
        Va = triple_gauge_vertex(eps3,  eps4,  mu5,  k1, k2, negq)
        Vc = triple_gauge_vertex(eps3c, eps4c, mu5c, k1, k2, negq)
        @test finish(_pair_trace(s_b, s_b, true)) == finish(tr * Va * Vc) * gs4
    end

    @testset "t-channel diagonal |M_t|² ≡ handbuilt D_tt" begin
        # Γ_t = γ^{:eps_4} (p̸1 − k̸1) γ^{:eps_3}  (quark exchange).
        p1mk1 = MomentumSum([(1//1, p1), (-1//1, k1)])
        gt  = DiracGamma[γ(eps4),  DiracGamma(MomSumSlot(p1mk1)), γ(eps3)]
        gtc = DiracGamma[γ(eps3c), DiracGamma(MomSumSlot(p1mk1)), γ(eps4c)]
        D_tt = dirac_trace(DiracGamma[GS(p1); gtc; GS(p2); gt])
        @test finish(_pair_trace(t_b, t_b, true)) == finish(D_tt) * gs4
    end

    @testset "u-channel diagonal |M_u|² ≡ handbuilt D_uu" begin
        # Γ_u = γ^{:eps_3} (p̸1 − k̸2) γ^{:eps_4}  (crossed quark exchange).
        p1mk2 = MomentumSum([(1//1, p1), (-1//1, k2)])
        gu  = DiracGamma[γ(eps3),  DiracGamma(MomSumSlot(p1mk2)), γ(eps4)]
        guc = DiracGamma[γ(eps4c), DiracGamma(MomSumSlot(p1mk2)), γ(eps3c)]
        D_uu = dirac_trace(DiracGamma[GS(p1); guc; GS(p2); gu])
        @test finish(_pair_trace(u_b, u_b, true)) == finish(D_uu) * gs4
    end

    @testset "cross-channel interference traces are non-trivial" begin
        # The 1×1 off-diagonal path (_line_trace) and the boson-factor
        # conjugation must both fire: s×t couples the off-line ggg factor
        # (conjugated) against a pure fermion line; t×u is fermion-only.
        @test !iszero(finish(_pair_trace(t_b, u_b, false)))
        @test !iszero(finish(_pair_trace(s_b, t_b, false)))
    end

    @testset "solve_tree_pipeline(qq̄→gg) assembles end-to-end" begin
        prob = CrossSectionProblem(qcd_model(m_q=:zero),
            [ExternalLeg(:q, p1, true, false), ExternalLeg(:q, p2, true, true)],
            [ExternalLeg(:g, k1, false, false), ExternalLeg(:g, k2, false, false)],
            10.0)
        r = solve_tree_pipeline(prob)
        @test r.n_emissions == 3
        @test r.amplitude_squared isa AlgSum
        @test !iszero(r.amplitude_squared)
        # NOTE: full qq̄→gg ≡ FeynCalc additionally needs colour algebra in
        # the pipeline (bead feynfeld-yewo) — see feynfeld-4xrh.
    end
end
