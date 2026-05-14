# Phase 18b-4 (bead feynfeld-m4o8): boson polarisation for external
# photons / gluons, validated on Compton scattering e γ → e γ.
#
# Compton is the clean infrastructure target: pure QED, one fermion
# line, two external photons, Feynman-gauge polarisation sum, no
# triple-gauge vertex, no colour. It exercises every piece m4o8 adds:
#   - external-boson ε^μ emission        (build_externals.pol_index)
#   - canonical :eps_<leg> relabelling   (emission_to_amplitude)
#   - 1-fermion-line bundle routing      (burnside_combine._pair_trace)
#   - ε–ε* polarisation sum plumbing     (combine_m_squared_burnside)
#
# Ground truth: Peskin & Schroeder Eq. (5.87), massless limit.
# Cited from refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/
# ElGa-ElGa.m, lines 103-105:
#   "|M̄|²/(2e⁴) = SP[p1,k2]/SP[p1,k1] + SP[p1,k1]/SP[p1,k2]
#                + 2m²(1/SP[p1,k1] - 1/SP[p1,k2])
#                + m⁴(1/SP[p1,k1] - 1/SP[p1,k2])²"
# with m = 0 (m_e=:zero pipeline): |M̄|²/(2e⁴) = p1·k2/p1·k1 + p1·k1/p1·k2.
#
# Two tree diagrams (s-channel + u-channel), both fermion_sign = +1, so
# they add — matching test/v2/test_compton.jl's handbuilt path.

using Test

using Feynfeld
using Feynfeld.QgrafPort: combine_m_squared_burnside, AmplitudeBundle,
                            emission_to_amplitude, compute_qg10_labels,
                            Partition, TopoState, MAX_V

@testset "Phase 18b-4: external-boson polarisation (Compton)" begin

    p1 = Momentum(:p1); p2 = Momentum(:p2)
    k1 = Momentum(:k1); k2 = Momentum(:k2)

    # ── Canonical :eps_<leg> indices the pipeline produces (ps1 fixes
    #    photon leg 2 = k1 → :eps_2, photon leg 4 = k2 → :eps_4). The
    #    conjugate amplitude relabels each index x → x_ (spin_sum.jl).
    eps_2  = LorentzIndex(:eps_2,  DimD()); eps_4  = LorentzIndex(:eps_4,  DimD())
    eps_2c = LorentzIndex(:eps_2_, DimD()); eps_4c = LorentzIndex(:eps_4_, DimD())
    γ(idx) = DiracGamma(LISlot(idx))

    # Per-bundle fermion-line gamma chains (middle of ū … u), pipeline
    # naming. s-channel: ū γ^{e4}(p̸1+k̸1)γ^{e2} u  (denom (p1+k1)²);
    # u-channel: ū γ^{e2}(p̸1-k̸2)γ^{e4} u          (denom (p1-k2)²).
    s_fwd = [(1//1, [γ(eps_4), GS(k1), γ(eps_2)]),
             (1//1, [γ(eps_4), GS(p1), γ(eps_2)])]
    u_fwd = [(1//1, [γ(eps_2), GS(p1), γ(eps_4)]),
             (-1//1,[γ(eps_2), GS(k2), γ(eps_4)])]
    # Conjugate chains: reverse + index relabel x → x_.
    s_cnj = [(1//1, [γ(eps_2c), GS(k1), γ(eps_4c)]),
             (1//1, [γ(eps_2c), GS(p1), γ(eps_4c)])]
    u_cnj = [(1//1, [γ(eps_4c), GS(p1), γ(eps_2c)]),
             (-1//1,[γ(eps_4c), GS(k2), γ(eps_2c)])]

    # T_ij = Σ_{fwd∈j} Σ_{cnj∈i} c_f c_c · Tr[p̸1 · Γ̄_cnj · p̸2 · Γ_fwd],
    # with the two Feynman-gauge photon polarisation sums applied and the
    # Lorentz indices contracted. Massless: completeness is just p̸.
    function Tblock(cnj_list, fwd_list)
        m = AlgSum()
        for (cf, gf) in fwd_list, (cc, gc) in cnj_list
            m = m + (cf * cc) * dirac_trace(DiracGamma[GS(p1); gc; GS(p2); gf])
        end
        m = m * polarization_sum(eps_2, eps_2c) * polarization_sum(eps_4, eps_4c)
        expand_scalar_product(contract(m))
    end
    T_ss = Tblock(s_cnj, s_fwd); T_su = Tblock(s_cnj, u_fwd)
    T_us = Tblock(u_cnj, s_fwd); T_uu = Tblock(u_cnj, u_fwd)
    handbuilt_trace = T_ss + T_su + T_us + T_uu

    prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
        [ExternalLeg(:e, p1, true, false), ExternalLeg(:gamma, k1, true, false)],
        [ExternalLeg(:e, p2, false, false), ExternalLeg(:gamma, k2, false, false)],
        10.0)

    @testset "Compton emits exactly 2 canonical orbits" begin
        result = solve_tree_pipeline(prob)
        @test result.n_emissions == 2
    end

    @testset "pipeline |M|² ≡ handbuilt trace-only (symbolic AlgSum ==)" begin
        # Acceptance criterion (1): qq̄→gg was descoped to Compton — see
        # bead feynfeld-m4o8 notes; full qq̄→gg → bead feynfeld-4xrh.
        result = solve_tree_pipeline(prob)
        @test result.amplitude_squared == handbuilt_trace
    end

    @testset "physics ground truth: |M̄|² == massless P&S Eq. (5.87)" begin
        # Re-validate the handbuilt reference IS Compton by restoring the
        # per-channel propagator denominators and the 1/4 spin+pol
        # average, then evaluating at a massless 2→2 point (s=4, t=-1,
        # u=-3). d_s = (p1+k1)² = s, d_u = (p1-k2)² = u.
        sv = 4//1; tv = -1//1; uv = -3//1
        ctx = sp_context(
            (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>0//1, (:k2,:k2)=>0//1,
            (:p1,:k1)=>sv//2,  (:p2,:k2)=>sv//2,
            (:p1,:k2)=>-uv//2, (:p2,:k1)=>-uv//2,
            (:p1,:p2)=>-tv//2, (:k1,:k2)=>-tv//2)
        ev(x) = (r = evaluate_sp(x; ctx); v = first(r.terms)[2];
                 v isa DimPoly ? evaluate_dim(v) : v)

        d_s = sv; d_u = uv
        m_sq = (1//4) * (ev(T_ss) // (d_s * d_s) + ev(T_su) // (d_s * d_u) +
                         ev(T_us) // (d_u * d_s) + ev(T_uu) // (d_u * d_u))

        # P&S 5.87, m → 0, e⁴ stripped: 2(p1·k2/p1·k1 + p1·k1/p1·k2).
        p1k1 = sv // 2; p1k2 = -uv // 2
        ps587_massless = 2 * (p1k2 // p1k1 + p1k1 // p1k2)
        @test m_sq == ps587_massless
    end

    @testset "combine_m_squared_burnside rejects mismatched boson_pols" begin
        # Two bundles of the same process must agree on their canonical
        # external-boson indices, else the polarisation sum is ambiguous.
        b1 = AmplitudeBundle(DiracExpr[], DiracExpr(alg(1)), AlgSum[], 1,
                             1//1, alg(1), [LorentzIndex(:eps_2, DimD())], alg(1))
        b2 = AmplitudeBundle(DiracExpr[], DiracExpr(alg(1)), AlgSum[], 1,
                             1//1, alg(1), [LorentzIndex(:eps_9, DimD())], alg(1))
        @test_throws ErrorException combine_m_squared_burnside([b1, b2], [1//1, 1//1])
    end

    @testset "Phase 18a/18b-1 regression: internal-boson processes unchanged" begin
        # ee→μμ: photon internal → empty boson_pols, no polarisation sum,
        # still a single orbit.
        prob_mumu = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
            [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)],
            [ExternalLeg(:mu, k1, false, false), ExternalLeg(:mu, k2, false, true)],
            10.0)
        r = solve_tree_pipeline(prob_mumu)
        @test r.n_emissions == 1
        @test !iszero(r.amplitude_squared)
    end
end
