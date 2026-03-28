# Spiral 2: Bhabha scattering e⁻(p1) + e⁺(p2) → e⁻(k1) + e⁺(k2)
#
# Ground truth: refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/ElAel-ElAel.m
# Lines 98-99:
# "knownResult = 2 SMP["e"]^4 (s^2+u^2)/t^2 + 4 SMP["e"]^4 u^2/(s t)
#              + 2 SMP["e"]^4 (t^2+u^2)/s^2"
# Verified CORRECT by FeynCalc (lines 100-102).
#
# Amplitude: M = M_t - M_s  (relative minus from identical fermion exchange)
#   M_t = (e²/t) [ū(k1) γ^α u(p1)] [v̄(p2) γ_α v(k2)]   (t-channel exchange)
#   M_s = (e²/s) [v̄(p2) γ^β u(p1)] [ū(k1) γ_β v(k2)]   (s-channel annihilation)
#
# Spin-averaged |M̄|²/e⁴ = (1/4) Σ_spins [|A_t/t - A_s/s|²]
#   Direct terms: product of two 4-gamma traces (separate fermion lines)
#   Interference: single 8-gamma trace (fermion lines reconnect)

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Spiral 2: Bhabha scattering" begin

    @testset "Bhabha |M̄|² massless vs FeynCalc" begin
        # Momenta
        p1 = Momentum(:p1)  # incoming electron
        p2 = Momentum(:p2)  # incoming positron
        k1 = Momentum(:k1)  # outgoing electron
        k2 = Momentum(:k2)  # outgoing positron

        # Massless Mandelstam: s + t + u = 0
        s_val = 10//1
        t_val = -3//1
        u_val = -(s_val + t_val)  # = -7

        # Scalar products (massless):
        # s = (p1+p2)² = 2p1·p2 = (k1+k2)² = 2k1·k2
        # t = (p1-k1)² = -2p1·k1 = (p2-k2)² = -2p2·k2
        # u = (p1-k2)² = -2p1·k2 = (p2-k1)² = -2p2·k1
        ctx = sp_context(
            (:p1,:p1) => 0//1, (:p2,:p2) => 0//1,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:p2) => s_val // 2,      # s/2 = 5
            (:p1,:k1) => -t_val // 2,     # -t/2 = 3/2
            (:p1,:k2) => -u_val // 2,     # -u/2 = 7/2
            (:p2,:k1) => -u_val // 2,     # -u/2 = 7/2
            (:p2,:k2) => -t_val // 2,     # -t/2 = 3/2
            (:k1,:k2) => s_val // 2,      # s/2 = 5
        )

        # ──── DIRECT TERM: Σ_spins |A_t|² ────
        # t-channel amplitude numerator: A_t = [ū(k1) γ^α u(p1)] [v̄(p2) γ_α v(k2)]
        # Conjugate: A_t* = [ū(p1) γ^{α'} u(k1)] [v̄(k2) γ_{α'} v(p2)]
        #
        # Spin sum line 1 (electron):
        #   [ū(p1) γ^{α'} u(k1)][ū(k1) γ^α u(p1)]
        #   → Tr[p̸₁ γ^{α'} k̸₁ γ^α]
        #
        # Spin sum line 2 (positron):
        #   [v̄(k2) γ_{α'} v(p2)][v̄(p2) γ_α v(k2)]
        #   → Tr[k̸₂ γ_{α'} p̸₂ γ_α]
        #
        # Ref: DESIGN.md anti-pattern #8: separate traces × multiply
        Tr_t1 = dirac_trace(DiracGamma[GS(p1), GAD(:alpha_), GS(k1), GAD(:alpha)])
        Tr_t2 = dirac_trace(DiracGamma[GS(k2), GAD(:alpha_), GS(p2), GAD(:alpha)])
        T_tt = Tr_t1 * Tr_t2

        # ──── DIRECT TERM: Σ_spins |A_s|² ────
        # s-channel amplitude numerator: A_s = [v̄(p2) γ^β u(p1)] [ū(k1) γ_β v(k2)]
        # Conjugate: A_s* = [ū(p1) γ^{β'} v(p2)] [v̄(k2) γ_{β'} u(k1)]
        #
        # Spin sum line 1:
        #   [ū(p1) γ^{β'} v(p2)][v̄(p2) γ^β u(p1)]
        #   → Tr[p̸₁ γ^{β'} p̸₂ γ^β]
        #
        # Spin sum line 2:
        #   [v̄(k2) γ_{β'} u(k1)][ū(k1) γ_β v(k2)]
        #   → Tr[k̸₂ γ_{β'} k̸₁ γ_β]
        Tr_s1 = dirac_trace(DiracGamma[GS(p1), GAD(:beta_), GS(p2), GAD(:beta)])
        Tr_s2 = dirac_trace(DiracGamma[GS(k2), GAD(:beta_), GS(k1), GAD(:beta)])
        T_ss = Tr_s1 * Tr_s2

        # ──── INTERFERENCE: Σ_spins A_s* A_t ────
        # When fermion lines reconnect across s and t channels, we get a SINGLE trace:
        #   Tr[p̸₁ γ^{β'} p̸₂ γ_α k̸₂ γ_{β'} k̸₁ γ^α]
        #
        # Derivation: M_s* × M_t spinor pairing forms one closed loop:
        #   ū(p1) → γ^{β'} → v(p2) → v̄(p2) → γ_α → v(k2)
        #   → v̄(k2) → γ_{β'} → u(k1) → ū(k1) → γ^α → u(p1)
        # After completeness relations (massless: u ū = p̸, v v̄ = p̸):
        #   Tr[p̸₁ γ^{β'} p̸₂ γ_α k̸₂ γ_{β'} k̸₁ γ^α]
        T_int = dirac_trace(DiracGamma[
            GS(p1), GAD(:beta_), GS(p2), GAD(:alpha),
            GS(k2), GAD(:beta_), GS(k1), GAD(:alpha)
        ])

        # By complex conjugation: Σ A_t* A_s = (Σ A_s* A_t)* = Σ A_s* A_t
        # (real in 4D with real momenta, confirmed by cyclic trace property)

        # ──── COMBINE: |M̄|²/e⁴ = (1/4)[T_tt/t² + T_ss/s² - 2·T_int/(st)] ────
        # Ref: ElAel-ElAel.m line 86: ExtraFactor -> 1/2^2 (spin average)
        # The minus sign: M = M_t - M_s, so cross term has -1 × (-1) for the two -M_s
        # Wait: |M_t - M_s|² = |M_t|² + |M_s|² - M_t*M_s - M_s*M_t
        # But T_int = Σ A_s* A_t, so the cross terms give -(1/(st))(T_int + T_int*)
        # = -2T_int/(st)
        total = 1//4 * (T_tt * (1 // t_val^2) + T_ss * (1 // s_val^2) -
                        2 * T_int * (1 // (s_val * t_val)))

        # ──── CONTRACT + EXPAND + EVALUATE ────
        contracted = contract(total; ctx)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx)

        # Extract numerical value
        @test length(result.terms) == 1  # single scalar term
        fk, coeff = first(result.terms)
        @test isempty(fk.factors)  # no remaining Lorentz factors

        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # ──── COMPARE TO FEYNCALC KNOWN RESULT ────
        # Ref: refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/ElAel-ElAel.m, lines 98-99
        # "|M̄|²/e⁴ = 2(s²+u²)/t² + 4u²/(st) + 2(t²+u²)/s²"
        # (coupling e⁴ stripped — we compute reduced matrix element)
        known = 2 * (s_val^2 + u_val^2) // t_val^2 +
                4 * u_val^2 // (s_val * t_val) +
                2 * (t_val^2 + u_val^2) // s_val^2

        @test pipeline_value == known
    end

    @testset "Bhabha |M̄|² second kinematic point" begin
        # Cross-check: different kinematics to catch accidental agreement.
        # s = 7, t = -2, u = -5
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        s2 = 7//1; t2 = -2//1; u2 = -5//1
        @assert s2 + t2 + u2 == 0

        ctx2 = sp_context(
            (:p1,:p1) => 0//1, (:p2,:p2) => 0//1,
            (:k1,:k1) => 0//1, (:k2,:k2) => 0//1,
            (:p1,:p2) => s2//2, (:p1,:k1) => -t2//2, (:p1,:k2) => -u2//2,
            (:p2,:k1) => -u2//2, (:p2,:k2) => -t2//2, (:k1,:k2) => s2//2,
        )

        # Same trace structures — only kinematics change
        Tr_t1 = dirac_trace(DiracGamma[GS(p1), GAD(:alpha_), GS(k1), GAD(:alpha)])
        Tr_t2 = dirac_trace(DiracGamma[GS(k2), GAD(:alpha_), GS(p2), GAD(:alpha)])
        Tr_s1 = dirac_trace(DiracGamma[GS(p1), GAD(:beta_), GS(p2), GAD(:beta)])
        Tr_s2 = dirac_trace(DiracGamma[GS(k2), GAD(:beta_), GS(k1), GAD(:beta)])
        T_int = dirac_trace(DiracGamma[
            GS(p1), GAD(:beta_), GS(p2), GAD(:alpha),
            GS(k2), GAD(:beta_), GS(k1), GAD(:alpha)])

        total = 1//4 * (Tr_t1 * Tr_t2 * (1 // t2^2) + Tr_s1 * Tr_s2 * (1 // s2^2) -
                        2 * T_int * (1 // (s2 * t2)))

        contracted = contract(total; ctx=ctx2)
        expanded = expand_scalar_product(contracted)
        result = evaluate_sp(expanded; ctx=ctx2)

        fk, coeff = first(result.terms)
        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # Ref: same formula, different numbers
        known2 = 2 * (s2^2 + u2^2) // t2^2 + 4 * u2^2 // (s2 * t2) +
                 2 * (t2^2 + u2^2) // s2^2

        @test pipeline_value == known2
    end
end
