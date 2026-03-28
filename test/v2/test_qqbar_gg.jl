# Spiral 3: qq̄ → gg (quark-antiquark annihilation to two gluons)
#
# Ground truth: refs/FeynCalc/FeynCalc/Examples/QCD/Tree/Mathematica/QQbar-GlGl.m
# Lines 103-104:
# "knownResult = (32/27) g_s^4 (t^2+u^2)/(t u) - (8/3) g_s^4 (t^2+u^2)/s^2"
# Verified CORRECT by FeynCalc vs Ellis, Stirling, Weber Table 7.1.
#
# Three diagrams:
#   t-channel: quark exchange, colour (T^b T^a)_{ji}
#   u-channel: crossed, colour (T^a T^b)_{ji}
#   s-channel: triple gluon vertex, colour T^c_{ji} f^{cab}
#
# Colour factors (N=3, summed over all external colours):
#   C_tt = C_uu = (N²-1)²/(4N) = 16/3
#   C_tu = C_ut = -(N²-1)/(4N) = -2/3
#   C_ss = N(N²-1)/2 = 12
#   C_ts, C_us = ±iN(N²-1)/4 = ±6i (pure imaginary → vanish at tree level)
#
# Spin+colour average: 1/(4 × N²) = 1/36

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Spiral 3: qq̄ → gg" begin

    @testset "qq̄→gg |M̄|² massless vs FeynCalc" begin
        # Momenta: q(p1) + q̄(p2) → g(k1) + g(k2)
        p1 = Momentum(:p1)
        p2 = Momentum(:p2)
        k1 = Momentum(:k1)
        k2 = Momentum(:k2)

        # Massless Mandelstam: s + t + u = 0
        s_val = 10//1
        t_val = -3//1
        u_val = -(s_val + t_val)  # = -7

        # Scalar products (massless)
        ctx = sp_context(
            (:p1,:p1) => 0//1, (:p2,:p2) => 0//1,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:p2) => s_val // 2,
            (:p1,:k1) => -t_val // 2,
            (:p1,:k2) => -u_val // 2,
            (:p2,:k1) => -u_val // 2,
            (:p2,:k2) => -t_val // 2,
            (:k1,:k2) => s_val // 2,
        )

        # ──── COLOUR FACTORS (analytical, N=3) ────
        # Ref: standard SU(N) Casimir identities
        C_tt = 16//3;  C_uu = 16//3
        C_tu = -2//3
        C_ss = 12//1
        # C_ts = -6i, C_us = +6i → vanish (Im × Re = 0 at tree level)

        # ──── DIRAC TRACES (spin + polarization summed) ────
        # After spin sum and pol sum (−g_{μμ'}, −g_{νν'}), each term is:
        #   D_XY = Tr[p̸₂ Γ_X p̸₁ Γ̃_Y] with contracted μ,ν indices
        #
        # t-channel Dirac: v̄(p2) γ^ν (p̸₁-k̸₁)/t γ^μ u(p1)
        # u-channel Dirac: v̄(p2) γ^μ (p̸₁-k̸₂)/u γ^ν u(p1)
        # s-channel Dirac: v̄(p2) γ^ρ u(p1) × (1/s) × V_{ρμν}

        # Build propagator MomentumSums
        p1mk1 = MomentumSum([(1//1, p1), (-1//1, k1)])  # p1 - k1
        p1mk2 = MomentumSum([(1//1, p1), (-1//1, k2)])  # p1 - k2

        # ──── t-channel and u-channel traces (same structure as Compton) ────
        # For |M_t|² spin+pol summed:
        # Tr[p̸₂ γ^ν (p̸₁-k̸₁) γ^μ p̸₁ γ_μ (p̸₁-k̸₁) γ_ν]
        # The pol sum contracts μ↔μ', ν↔ν'. We relabel conjugate indices → same.
        # Then contract() handles repeated indices.

        # Build the trace for each diagram pair using the Compton pattern.
        # Trace chain for diagram X with conjugate Y:
        # Tr[p̸₂ Γ_X p̸₁ Γ̃_Y] where Γ̃_Y is reversed + relabeled

        # t-channel gamma chain (between completeness insertions):
        # Γ_t = γ^ν (p̸₁-k̸₁) γ^μ (the /t denominator handled separately)
        gamma_t = DiracGamma[GAD(:nu), DiracGamma(MomSumSlot(p1mk1)), GAD(:mu)]
        gamma_t_conj = DiracGamma[GAD(:mu), DiracGamma(MomSumSlot(p1mk1)), GAD(:nu)]

        # u-channel: γ^μ (p̸₁-k̸₂) γ^ν
        gamma_u = DiracGamma[GAD(:mu), DiracGamma(MomSumSlot(p1mk2)), GAD(:nu)]
        gamma_u_conj = DiracGamma[GAD(:nu), DiracGamma(MomSumSlot(p1mk2)), GAD(:mu)]

        # Helper: trace Tr[p̸₂ Γ_X p̸₁ Γ̃_Y] for massless fermions
        function dirac_trace_XY(gamma_X, gamma_Y_conj)
            gs = DiracGamma[GS(p2); gamma_X; GS(p1); gamma_Y_conj]
            dirac_trace(gs)
        end

        # Use PRIMED indices for conjugate, apply PHYSICAL pol sum later
        # (QCD requires physical pol sums — Ward identity only holds for full amplitude)
        gamma_t_conj_p = DiracGamma[GAD(:mu_), DiracGamma(MomSumSlot(p1mk1)), GAD(:nu_)]
        gamma_u_conj_p = DiracGamma[GAD(:nu_), DiracGamma(MomSumSlot(p1mk2)), GAD(:mu_)]

        D_tt_raw = dirac_trace_XY(gamma_t, gamma_t_conj_p)
        D_uu_raw = dirac_trace_XY(gamma_u, gamma_u_conj_p)
        D_tu_raw = dirac_trace_XY(gamma_t, gamma_u_conj_p)

        # ──── s-channel ────
        # Triple gluon vertex V_{ρμν}(-q,k1,k2) with q=k1+k2 (all outgoing):
        # V = g_{ρμ}(-2k1-k2)_ν + g_{μν}(k1-k2)_ρ + g_{νρ}(k1+2k2)_μ
        rho = LorentzIndex(:rho, DimD())
        sig = LorentzIndex(:sig, DimD())
        mu_idx = LorentzIndex(:mu, DimD())
        nu_idx = LorentzIndex(:nu, DimD())
        mu_p = LorentzIndex(:mu_, DimD())
        nu_p = LorentzIndex(:nu_, DimD())

        m2k1mk2 = MomentumSum([(-2//1, k1), (-1//1, k2)])
        k1mk2   = MomentumSum([(1//1, k1), (-1//1, k2)])
        k1p2k2  = MomentumSum([(1//1, k1), (2//1, k2)])

        function vertex_term(g1, g2, ms, mi)
            result = AlgSum()
            for (c, m) in ms.terms
                result = result + c * alg(pair(g1, g2)) * alg(pair(mi, m))
            end
            result
        end

        V_amp = vertex_term(rho, mu_idx, m2k1mk2, nu_idx) +
                vertex_term(mu_idx, nu_idx, k1mk2, rho) +
                vertex_term(nu_idx, rho, k1p2k2, mu_idx)

        V_conj = vertex_term(sig, mu_p, m2k1mk2, nu_p) +
                 vertex_term(mu_p, nu_p, k1mk2, sig) +
                 vertex_term(nu_p, sig, k1p2k2, mu_p)

        qt_ss = dirac_trace(DiracGamma[GS(p2), GAD(:rho), GS(p1), GAD(:sig)])
        D_ss_raw = qt_ss * V_amp * V_conj

        # ──── COMBINE ALL with colour factors BEFORE pol sum ────
        total_raw = (1//36) * (
            (C_tt // t_val^2) * D_tt_raw +
            (C_uu // u_val^2) * D_uu_raw +
            (C_tu // (t_val * u_val)) * D_tu_raw +
            (C_tu // (t_val * u_val)) * dirac_trace_XY(gamma_u, gamma_t_conj_p) +
            (C_ss // s_val^2) * D_ss_raw
        )

        # ──── PHYSICAL POLARIZATION SUMS (axial gauge) ────
        # Ref: FeynCalc DoPolarizationSums[#,k1,k2] convention
        # For k1: Σ ε^μ ε^{μ'*} = -g^{μμ'} + (k1^μ k2^{μ'} + k2^μ k1^{μ'})/(k1·k2)
        # For k2: Σ ε^ν ε^{ν'*} = -g^{νν'} + (k2^ν k1^{ν'} + k1^ν k2^{ν'})/(k1·k2)
        pol1 = polarization_sum(mu_idx, mu_p, k1, k2; ctx)
        pol2 = polarization_sum(nu_idx, nu_p, k2, k1; ctx)

        # Multiply total by both pol sums, then contract
        total_pol = total_raw * pol1 * pol2

        contracted = contract(total_pol; ctx)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx)

        # Extract numerical value
        @test length(result.terms) == 1
        fk, coeff = first(result.terms)
        @test isempty(fk.factors)

        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # ──── COMPARE TO FEYNCALC KNOWN RESULT ────
        # Ref: refs/FeynCalc/FeynCalc/Examples/QCD/Tree/Mathematica/QQbar-GlGl.m, lines 103-104
        # "|M̄|²/g_s⁴ = (32/27)(t²+u²)/(tu) - (8/3)(t²+u²)/s²"
        known = (32//27) * (t_val^2 + u_val^2) // (t_val * u_val) -
                (8//3) * (t_val^2 + u_val^2) // s_val^2

        @test pipeline_value == known
    end
end
