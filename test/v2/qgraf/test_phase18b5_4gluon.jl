# Phase 18b-5 (bead feynfeld-awtt): 4-gluon contact vertex Lorentz factor
# wired into the qgraf amplitude pipeline. Extends asp9's all-boson
# vertex dispatch in `_vertex_factor_at` from vdeg=3 to vdeg=4 by
# introducing a sibling `quadruple_gauge_vertex` (Layer 4) and a
# `_quadruple_boson_vertex_factor` helper in `build_vertices`.
#
# Ground truth: Peskin & Schroeder Eq. (16.5)-(16.6), four-gluon
# Feynman rule. Independent cross-check: refs/FeynCalc/.../Feynman/
# GluonVertex.m:108-118 — three colour-paired Lorentz tensors,
#
#   -i g_s² [ f^{abe}f^{cde} (g_{μλ}g_{νσ} − g_{μσ}g_{νλ})
#           + f^{ace}f^{bde} (g_{μν}g_{λσ} − g_{μσ}g_{νλ})
#           + f^{ade}f^{bce} (g_{μν}g_{λσ} − g_{μλ}g_{νσ}) ].
#
# Colour structures (`f^{abe}f^{cde}`, etc.) are deferred to bead
# `feynfeld-yewo` (no colour field in `AmplitudeBundle` yet). Like
# `triple_gauge_vertex` (asp9), the Layer-4 helper here returns the
# bare sum of three Lorentz tensors; the colour pairings are restored
# when yewo plumbs colour through the pipeline. Per the orchestrator
# brief, that's an acceptable systematic deferral parallel to ggg's.
#
# Acceptance criteria (from bead feynfeld-awtt):
#   (1) gg→gg tree pipeline emission count = 4.
#   (2) symbolic equality of pipeline `AmplitudeBundle` factors against
#       handbuilt references built from the cited `quadruple_gauge_vertex`
#       (contact) and `triple_gauge_vertex` (s/t/u).
#   (3) regression: qq̄→gg (asp9) still green — covered by
#       `test_triple_gauge_vertex.jl`.

using Test

using Feynfeld
using Feynfeld.QgrafPort: AmplitudeBundle, emission_to_amplitude,
                            enumerate_topology_automorphisms,
                            is_emission_canonical, _foreach_emission

@testset "Phase 18b-5: 4-gluon vertex (gg→gg through the pipeline)" begin

    p1 = Momentum(:p1); p2 = Momentum(:p2)
    k1 = Momentum(:k1); k2 = Momentum(:k2)

    # ── Layer-4 unit: `quadruple_gauge_vertex` is a pure algebraic helper.
    @testset "quadruple_gauge_vertex is the 3-term Lorentz tensor" begin
        m1 = LorentzIndex(:m1, DimD()); m2 = LorentzIndex(:m2, DimD())
        m3 = LorentzIndex(:m3, DimD()); m4 = LorentzIndex(:m4, DimD())
        # Verbatim P&S 16.5-16.6 / FeynCalc GluonVertex.m:116-118
        # Lorentz tensors (colour stripped, the three terms are summed
        # because colour algebra is the yewo deferral):
        #   T12_34 = g_{m1,m3}g_{m2,m4} − g_{m1,m4}g_{m2,m3}
        #   T13_24 = g_{m1,m2}g_{m3,m4} − g_{m1,m4}g_{m2,m3}
        #   T14_23 = g_{m1,m2}g_{m3,m4} − g_{m1,m3}g_{m2,m4}
        T12_34 = alg(pair(m1, m3)) * alg(pair(m2, m4)) -
                 alg(pair(m1, m4)) * alg(pair(m2, m3))
        T13_24 = alg(pair(m1, m2)) * alg(pair(m3, m4)) -
                 alg(pair(m1, m4)) * alg(pair(m2, m3))
        T14_23 = alg(pair(m1, m2)) * alg(pair(m3, m4)) -
                 alg(pair(m1, m3)) * alg(pair(m2, m4))
        @test quadruple_gauge_vertex(m1, m2, m3, m4) == T12_34 + T13_24 + T14_23
        # Symmetric under simultaneous swap of two legs: swapping (1↔2) and
        # (3↔4) leaves the tensor invariant (each pair just renames into
        # itself; T12_34 stays, T13_24 ↔ T14_23 swap, sum unchanged).
        @test quadruple_gauge_vertex(m1, m2, m3, m4) ==
              quadruple_gauge_vertex(m2, m1, m4, m3)
    end

    # ── Pipeline: gg→gg emits 4 canonical orbits (contact + s + t + u).
    function gggg_bundles()
        model = qcd_model(m_q=:zero)
        legs = [ExternalLeg(:g, p1, true,  false),
                ExternalLeg(:g, p2, true,  false),
                ExternalLeg(:g, k1, false, false),
                ExternalLeg(:g, k2, false, false)]
        ext_exp = Feynfeld._qgraf_ext_labels(model, legs)
        pm = [l.momentum for l in legs]
        pa = Bool[l.antiparticle for l in legs]
        bundles = AmplitudeBundle[]
        _foreach_emission(model, [:g, :g], [:g, :g];
                          loops=0, ext_exp=ext_exp) do st, lb, ps1, pmap
            is_emission_canonical(st, lb, enumerate_topology_automorphisms(st),
                                  ps1, pmap) || return
            push!(bundles, emission_to_amplitude(st, lb, ps1, pmap, model;
                    physical_moms=pm, n_inco=2, phys_anti=pa))
        end
        bundles
    end

    @testset "gg→gg emits exactly 4 canonical orbits" begin
        bundles = gggg_bundles()
        @test length(bundles) == 4
    end

    # ── The contact (vdeg=4) bundle is the one with zero propagators.
    # Cross-check: every cubic (vdeg=3) bundle has exactly one internal
    # gluon → exactly one denominator. The contact bundle has none.
    @testset "exactly one bundle has zero denoms (the contact diagram)" begin
        bundles = gggg_bundles()
        @test count(b -> isempty(b.denoms), bundles) == 1
        @test count(b -> length(b.denoms) == 1, bundles) == 3
        for b in bundles
            @test isempty(b.line_chains)     # no fermion lines for gg→gg
            @test b.boson_pols == [LorentzIndex(:eps_1, DimD()),
                                    LorentzIndex(:eps_2, DimD()),
                                    LorentzIndex(:eps_3, DimD()),
                                    LorentzIndex(:eps_4, DimD())]
        end
    end

    # ── Contact bundle: boson_factor ≡ quadruple_gauge_vertex on the
    # external pol indices. All four legs carry index `:eps_<phys_leg>`
    # at this stage (emission_to_amplitude canonicalises `:mu_l_<i>` for
    # external bosons). Identity ps1 here: leg slot i = physical leg i,
    # so `:eps_i` for i ∈ {1,2,3,4}.
    @testset "contact bundle: boson_factor ≡ quadruple_gauge_vertex" begin
        bundles = gggg_bundles()
        contact = bundles[findfirst(b -> isempty(b.denoms), bundles)]
        eps_(i) = LorentzIndex(Symbol(:eps_, i), DimD())
        @test contact.boson_factor ==
              quadruple_gauge_vertex(eps_(1), eps_(2), eps_(3), eps_(4))
        # No fermion line → amplitude is the boson_factor lifted to DiracExpr.
        @test contact.amplitude == DiracExpr(contact.boson_factor)
    end

    # ── solve_tree_pipeline end-to-end: 4 emissions, AlgSum result. As
    # with asp9 (qq̄→gg), the trace-only |M|² is NOT yet comparable to
    # FeynCalc per-channel because colour algebra is still deferred
    # (bead feynfeld-yewo → feynfeld-4xrh). The assertion guards the
    # pipeline assembly does not regress on the all-boson 4-vertex path.
    @testset "solve_tree_pipeline(gg→gg) assembles end-to-end" begin
        prob = CrossSectionProblem(qcd_model(m_q=:zero),
            [ExternalLeg(:g, p1, true,  false), ExternalLeg(:g, p2, true,  false)],
            [ExternalLeg(:g, k1, false, false), ExternalLeg(:g, k2, false, false)],
            10.0)
        r = solve_tree_pipeline(prob)
        @test r.n_emissions == 4
        @test r.amplitude_squared isa AlgSum
        @test !iszero(r.amplitude_squared)
    end
end
