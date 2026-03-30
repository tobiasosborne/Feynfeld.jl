# Pipeline tests: all processes through Model -> Rules -> Channels -> Amplitude -> Algebra.
#
# Validates that pipeline-generated amplitudes match existing hand-built results.
# Each test uses tree_channels + build_amplitude instead of direct chain construction.

using Test
@isdefined(FeynfeldX) || include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Pipeline: all processes" begin

    # ======== Bhabha: e-(p1) e+(p2) -> e-(k1) e+(k2) ========
    @testset "Bhabha e-e+ -> e-e+ via pipeline" begin
        model = qed_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
        outgoing = [ExternalLeg(:e, k1, false, false), ExternalLeg(:e, k2, false, true)]

        channels = tree_channels(model, rules, incoming, outgoing)
        @test length(channels) == 2

        # Identify s and t channels
        ch_s = first(c for c in channels if c.channel == :s)
        ch_t = first(c for c in channels if c.channel == :t)

        amp_s = build_amplitude(ch_s, rules, model)
        amp_t = build_amplitude(ch_t, rules, model)

        # Direct terms: Tr x Tr for each channel
        T_ss = spin_sum_amplitude_squared(amp_s[1], amp_s[2])
        T_tt = spin_sum_amplitude_squared(amp_t[1], amp_t[2])

        # Interference: reconnected single trace
        T_int = spin_sum_interference(amp_s, amp_t)

        # Kinematics: s=10, t=-3, u=-7
        s_val = 10//1; t_val = -3//1; u_val = -(s_val + t_val)

        ctx = sp_context(
            (:p1,:p1) => 0//1, (:p2,:p2) => 0//1,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:p2) => s_val // 2, (:p1,:k1) => -t_val // 2,
            (:p1,:k2) => -u_val // 2, (:p2,:k1) => -u_val // 2,
            (:p2,:k2) => -t_val // 2, (:k1,:k2) => s_val // 2,
        )

        # Combine: |M-bar|^2/e^4 = (1/4)[T_ss/s^2 + T_tt/t^2 - 2*T_int/(s*t)]
        # Ref: M = M_t - M_s, so cross term is -2*Re(M_s* M_t)
        total = 1//4 * (T_ss * (1 // s_val^2) + T_tt * (1 // t_val^2) -
                        2 * T_int * (1 // (s_val * t_val)))

        contracted = contract(total; ctx)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx)
        fk, coeff = first(result.terms)
        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # Ref: refs/FeynCalc/.../ElAel-ElAel.m, lines 98-99
        # "|M-bar|^2/e^4 = 2(s^2+u^2)/t^2 + 4u^2/(st) + 2(t^2+u^2)/s^2"
        known = 2 * (s_val^2 + u_val^2) // t_val^2 +
                4 * u_val^2 // (s_val * t_val) +
                2 * (t_val^2 + u_val^2) // s_val^2

        @test pipeline_value == known
    end

    @testset "Bhabha second kinematic point" begin
        model = qed_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
        outgoing = [ExternalLeg(:e, k1, false, false), ExternalLeg(:e, k2, false, true)]

        channels = tree_channels(model, rules, incoming, outgoing)
        ch_s = first(c for c in channels if c.channel == :s)
        ch_t = first(c for c in channels if c.channel == :t)
        amp_s = build_amplitude(ch_s, rules, model)
        amp_t = build_amplitude(ch_t, rules, model)

        T_ss = spin_sum_amplitude_squared(amp_s[1], amp_s[2])
        T_tt = spin_sum_amplitude_squared(amp_t[1], amp_t[2])
        T_int = spin_sum_interference(amp_s, amp_t)

        s2 = 7//1; t2 = -2//1; u2 = -5//1
        ctx2 = sp_context(
            (:p1,:p1) => 0//1, (:p2,:p2) => 0//1,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:p2) => s2//2, (:p1,:k1) => -t2//2, (:p1,:k2) => -u2//2,
            (:p2,:k1) => -u2//2, (:p2,:k2) => -t2//2, (:k1,:k2) => s2//2,
        )

        total = 1//4 * (T_ss * (1//s2^2) + T_tt * (1//t2^2) - 2*T_int*(1//(s2*t2)))
        contracted = contract(total; ctx=ctx2)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx=ctx2)
        fk, coeff = first(result.terms)
        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        known2 = 2*(s2^2+u2^2)//t2^2 + 4*u2^2//(s2*t2) + 2*(t2^2+u2^2)//s2^2
        @test pipeline_value == known2
    end

    # ======== Compton: e-(p1) gamma(k1) -> e-(p2) gamma(k2) ========
    @testset "Compton via pipeline" begin
        model = qed_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        # Massive electrons (m = 1, m^2 = 1)
        incoming = [ExternalLeg(:e, p1, true, false, 1//1), ExternalLeg(:gamma, k1, true, false)]
        outgoing = [ExternalLeg(:e, p2, false, false, 1//1), ExternalLeg(:gamma, k2, false, false)]

        channels = tree_channels(model, rules, incoming, outgoing)
        @test length(channels) == 2  # s-channel + u-channel (fermion exchange)

        # Build amplitude for each channel: (chain_mom, chain_mass, mass)
        amps = [build_amplitude(ch, rules, model) for ch in channels]

        # Massive Compton: m^2 = 1, s = 5, u = -2, t = 2m^2 - s - u = -1
        m2_val = 1//1; m_val = 1//1
        s_val = 5//1; u_val = -2//1
        t_val = 2 * m2_val - s_val - u_val

        p1k1 = (s_val - m2_val) // 2   # = 2
        p1k2 = (m2_val - u_val) // 2   # = 3/2

        ctx = sp_context(
            (:p1,:p1) => m2_val, (:p2,:p2) => m2_val,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:k1) => p1k1, (:p1,:k2) => p1k2,
            (:p2,:k1) => p1k2, (:p2,:k2) => p1k1,
            (:p1,:p2) => (2*m2_val - t_val)//2,
            (:k1,:k2) => -t_val//2,
        )

        # Propagator denominators
        denoms = [ch.channel == :s ? s_val - m2_val : u_val - m2_val for ch in channels]

        # Helper: compute Tr for a pair of diagram components
        # chain_fwd provides Gamma_j, chain_conj provides Gamma_i (conjugated)
        function trace_pair(chain_fwd, chain_conj)
            _cross_trace = FeynfeldX._cross_line_trace(chain_fwd, chain_conj)
        end

        # For each (i,j) diagram pair, compute the full trace
        # including propagator mass expansion:
        # T_ij = Tr[(p2+m)Gamma_j_conj(p1+m)Gamma_i]
        # where Gamma_i = gamma_i_mom + m*gamma_i_mass
        # Expanding: 4 sub-traces weighted by mass powers
        function full_trace(amp_i, amp_j)
            ch_i_mom, ch_i_mass, mi = amp_i
            ch_j_mom, ch_j_mass, mj = amp_j
            result = trace_pair(ch_i_mom, ch_j_mom)
            result = result + mj * trace_pair(ch_i_mom, ch_j_mass)
            result = result + mi * trace_pair(ch_i_mass, ch_j_mom)
            result = result + mi * mj * trace_pair(ch_i_mass, ch_j_mass)
            result
        end

        # Compute all T_{ij}, apply pol sums (relabel primed -> unprimed), contract
        total = AlgSum()
        for (i, amp_i) in enumerate(amps)
            for (j, amp_j) in enumerate(amps)
                T = full_trace(amp_i, amp_j)
                total = total + (1 // (denoms[i] * denoms[j])) * T
            end
        end
        total = 1//4 * total  # spin+pol average: 1/(2 spins × 2 pols)

        # Polarization sums: relabel primed indices -> unprimed, then contract.
        # Indices are DimD (from build_amplitude).
        mu1  = LorentzIndex(:mu_k1, DimD())
        mu1_ = LorentzIndex(:mu_k1_, DimD())
        mu2  = LorentzIndex(:mu_k2, DimD())
        mu2_ = LorentzIndex(:mu_k2_, DimD())
        pol_summed = substitute_index(total, mu1_, mu1)
        pol_summed = substitute_index(pol_summed, mu2_, mu2)

        contracted = contract(pol_summed; ctx)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx)

        @test length(result.terms) == 1
        fk, coeff = first(result.terms)
        @test isempty(fk.factors)
        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # Ref: P&S Eq. (5.87), FeynCalc ElGa-ElGa.m lines 103-105
        # |M-bar|^2/(2e^4) = p1.k2/p1.k1 + p1.k1/p1.k2
        #   + 2m^2(1/p1.k1 - 1/p1.k2) + m^4(1/p1.k1 - 1/p1.k2)^2
        ps587 = 2 * (p1k2 // p1k1 + p1k1 // p1k2 +
                     2 * m2_val * (1 // p1k1 - 1 // p1k2) +
                     m2_val^2 * (1 // p1k1 - 1 // p1k2)^2)

        @test pipeline_value == ps587
    end

    # ======== qq-bar -> gg: 3 channels (t, u fermion + s gluon) ========
    @testset "qq-bar -> gg via pipeline" begin
        model = qcd_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:q, p1, true, false), ExternalLeg(:q, p2, true, true)]
        outgoing = [ExternalLeg(:g, k1, false, false), ExternalLeg(:g, k2, false, false)]

        channels = tree_channels(model, rules, incoming, outgoing)
        @test length(channels) == 3  # s (gluon exchange), t, u (quark exchange)

        s_val = 10//1; t_val = -3//1; u_val = -(s_val + t_val)
        ctx = sp_context(
            (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>0//1, (:k2,:k2)=>0//1,
            (:p1,:p2)=>s_val//2, (:p1,:k1)=>-t_val//2, (:p1,:k2)=>-u_val//2,
            (:p2,:k1)=>-u_val//2, (:p2,:k2)=>-t_val//2, (:k1,:k2)=>s_val//2)

        # Physical (axial gauge) polarization sums with reference momenta k1<->k2.
        # Index names match build_amplitude's convention: :mu_k1, :mu_k2.
        mi = LorentzIndex(:mu_k1, DimD()); mp = LorentzIndex(:mu_k1_, DimD())
        ni = LorentzIndex(:mu_k2, DimD()); np = LorentzIndex(:mu_k2_, DimD())
        P1 = polarization_sum(mi, mp, k1, k2; ctx)
        P2 = polarization_sum(ni, np, k2, k1; ctx)

        function eval_D(expr)
            c = contract(expr * P1 * P2; ctx)
            e = expand_scalar_product(c)
            r = evaluate_sp(e; ctx)
            v = first(r.terms)[2]
            v isa DimPoly ? evaluate_dim(v) : v
        end

        # ---- t/u channels: fermion exchange via pipeline ----
        ch_t = first(c for c in channels if c.channel == :t)
        ch_u = first(c for c in channels if c.channel == :u)

        amp_t = build_amplitude(ch_t, rules, model)
        amp_u = build_amplitude(ch_u, rules, model)

        # For massless quarks, only chain_mom contributes
        function fermion_trace(amp_fwd, amp_conj)
            FeynfeldX._cross_line_trace(amp_fwd[1], amp_conj[1])
        end

        D_tt = eval_D(fermion_trace(amp_t, amp_t))
        D_uu = eval_D(fermion_trace(amp_u, amp_u))
        D_tu = eval_D(fermion_trace(amp_t, amp_u))
        D_ut = eval_D(fermion_trace(amp_u, amp_t))

        # ---- s-channel: via pipeline (gauge exchange) ----
        ch_s = first(c for c in channels if c.channel == :s)
        chain_s, vtx_s = build_amplitude(ch_s, rules, model)

        # Spin sum: Tr[p̸_out γ^ρ p̸_in γ^σ] from the fermion chain
        rho_s = LorentzIndex(:rho_s, DimD()); sig_s = LorentzIndex(:sig_s, DimD())
        qt_ss = dirac_trace(DiracGamma[GS(p2), DiracGamma(LISlot(rho_s)),
                                        GS(p1), DiracGamma(LISlot(sig_s))])

        # Conjugate vertex: same structure with sig_s and primed external indices
        neg_q = MomentumSum([(-1//1, p1), (-1//1, p2)])
        Vc = triple_gauge_vertex(sig_s, mp, np, neg_q, k1, k2)
        D_ss = eval_D(qt_ss * vtx_s * Vc)

        # ---- Combine with SU(3) colour factors ----
        # Ref: refs/FeynCalc/.../QQbar-GlGl.m, lines 98-105
        C_tt = 16//3; C_uu = 16//3; C_tu = -2//3; C_ss = 12//1

        total = (1//36) * (
            C_tt * D_tt // t_val^2 +
            C_uu * D_uu // u_val^2 +
            C_tu * D_tu // (t_val * u_val) +
            C_tu * D_ut // (t_val * u_val) +
            C_ss * (-D_ss) // s_val^2)  # -D_ss from relative phase

        # Ref: "|M-bar|^2/g_s^4 = (32/27)(t^2+u^2)/(tu) - (8/3)(t^2+u^2)/s^2"
        known = (32//27) * (t_val^2 + u_val^2) // (t_val * u_val) -
                (8//3) * (t_val^2 + u_val^2) // s_val^2

        @test total == known
    end

    # ======== EW: e+e- -> W+W- pipeline (amplitude building) ========
    @testset "ee -> WW via pipeline" begin
        model = ew_model()
        rules = feynman_rules(model)
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
        outgoing = [ExternalLeg(:W, k1, false, false), ExternalLeg(:W, k2, false, true)]

        channels = tree_channels(model, rules, incoming, outgoing)

        # Channel enumeration
        @test length(channels) == 4  # s-γ, s-Z, t-ν, u-ν
        s_channels = [c for c in channels if c.channel == :s]
        t_channels = [c for c in channels if c.channel == :t]
        @test length(s_channels) == 2
        @test length(t_channels) == 1

        # Build ALL channel amplitudes through the pipeline
        ch_sg = first(c for c in channels if c.channel == :s && c.exchanged == :gamma)
        ch_sZ = first(c for c in channels if c.channel == :s && c.exchanged == :Z)
        ch_t  = first(c for c in channels if c.channel == :t)

        amp_sg = build_amplitude(ch_sg, rules, model)
        amp_sZ = build_amplitude(ch_sZ, rules, model)
        amp_t  = build_amplitude(ch_t, rules, model)

        # s-γ: gauge exchange → (DiracExpr fermion chain, AlgSum gauge vertex)
        @test amp_sg isa Tuple{DiracExpr, AlgSum}
        @test length(amp_sg[1].terms) == 1   # QED vertex: single γ^μ term

        # s-Z: gauge exchange → chiral vertex has 2 terms (g_V γ^μ - g_A γ5 γ^μ)
        @test amp_sZ isa Tuple{DiracExpr, AlgSum}
        @test length(amp_sZ[1].terms) == 2   # chiral: 2 terms

        # t-ν: fermion exchange → chain with GA7 chiral projectors
        @test amp_t isa Tuple{DiracChain, DiracChain, Rational{Int}}
        # Verify chiral projectors present: GA7() at both vertices
        t_elems = amp_t[1].elements
        @test any(e -> e isa DiracGamma{ProjMSlot}, t_elems)  # (1-γ5)/2

        # ---- Evaluate s-γ channel at a specific kinematic point ----
        # Use clean Rationals: s=8 (above threshold 4), t=-1, u=4-8-(-1)=-3
        # In M_W units: k1²=k2²=1 (M_W²), p1²=p2²=0 (massless electrons)
        s_val = 8//1; t_val = -1//1; u_val = 4//1 - s_val - t_val  # = -3

        ctx = sp_context(
            (:p1,:p1)=>0//1, (:p2,:p2)=>0//1,
            (:k1,:k1)=>1//1, (:k2,:k2)=>1//1,
            (:p1,:p2)=>s_val//2, (:p1,:k1)=>-t_val//2,
            (:p1,:k2)=>-u_val//2, (:p2,:k1)=>-u_val//2,
            (:p2,:k2)=>-t_val//2, (:k1,:k2)=>(s_val-2)//2)

        # Spin-summed fermion trace for s-γ
        chain_sg, vtx_sg = amp_sg
        tr_sg = spin_sum_amplitude_squared(chain_sg, chain_sg)  # single line → trace²? no...
        # Actually: gauge exchange has one fermion line → single trace
        tr_sg = FeynfeldX._single_line_trace(chain_sg)

        # Contract gauge vertex squared with W polarization sums
        rho_s = LorentzIndex(:rho_s, DimD())
        sig_s = LorentzIndex(:sig_s, DimD())
        neg_q = MomentumSum([(-1//1, p1), (-1//1, p2)])
        vtx_conj = triple_gauge_vertex(sig_s,
            LorentzIndex(:mu_k1_, DimD()), LorentzIndex(:mu_k2_, DimD()),
            neg_q, k1, k2)

        # W polarization sums (massive)
        mi  = LorentzIndex(:mu_k1, DimD());  mp = LorentzIndex(:mu_k1_, DimD())
        ni  = LorentzIndex(:mu_k2, DimD());  np = LorentzIndex(:mu_k2_, DimD())
        P1 = polarization_sum_massive(mi, mp, k1, 1//1)   # M_W²=1 in M_W units
        P2 = polarization_sum_massive(ni, np, k2, 1//1)

        # Full contraction
        full = tr_sg * vtx_sg * vtx_conj * P1 * P2
        c = contract(full; ctx)
        e = expand_scalar_product(c)
        r = evaluate_sp(e; ctx)
        val = first(r.terms)[2]
        result = val isa DimPoly ? Float64(evaluate_dim(val)) : Float64(val)

        @test isfinite(result)
        @test result != 0.0  # nontrivial result from pipeline
    end
end
