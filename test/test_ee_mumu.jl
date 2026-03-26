# End-to-end tracer bullet: e+e- → μ+μ- at tree level
#
# Validates Peskin & Schroeder Eq (5.10):
#   dσ/dΩ = α²/(4s) (1 + cos²θ)     (CM frame, massless limit)
#
# Equivalent to:
#   (1/4) Σ_spins |M|² = 2e⁴(t² + u²)/s²
#
# where s,t,u are Mandelstam variables with s+t+u=0 (massless).
#
# Process: e⁻(p) e⁺(p') → μ⁻(k) μ⁺(k'), s-channel photon exchange
#
# Source: FeynCalc Examples/QED/Tree/Mathematica/ElAel-MuAmu.m

using Feynfeld: fermion_spin_sum, contract, evaluate_sp, expand_scalar_product,
    AlgSum, AlgTerm, alg, alg_zero, alg_scalar, is_scalar,
    Spinor, DiracChain, dot, GA, GS,
    Pair, SP, FV, MT, LorentzIndex, Momentum,
    SPContext, set_sp, dirac_trace_alg

@testset "e+e- → μ+μ- tree-level (P&S 5.10)" begin

    # ── Step 1: Set up amplitude structure ─────────────────────────
    # M ∝ [ū(k) γ^μ v(k')] × (-g_{μν}/s) × [v̄(p') γ^ν u(p)]
    # Use same indices μ,ν in both lines (propagator already contracted)

    k  = Momentum(:k)     # μ⁻
    kp = Momentum(:kp)    # μ⁺
    p  = Momentum(:p)     # e⁻
    pp = Momentum(:pp)    # e⁺

    # ── Step 2: Spin sums via FermionSpinSum ───────────────────────

    # Muon line: [ū(k) γ^μ v(k')] [v̄(k') γ^ν u(k)]
    mu_M  = dot(Spinor(:ubar, k), GA(:μ), Spinor(:v, kp))
    mu_Mc = dot(Spinor(:vbar, kp), GA(:ν), Spinor(:u, k))
    tr_mu = fermion_spin_sum(mu_M, mu_Mc)

    # Electron line: [v̄(p') γ^μ u(p)] [ū(p) γ^ν v(p')]
    e_M  = dot(Spinor(:vbar, pp), GA(:μ), Spinor(:u, p))
    e_Mc = dot(Spinor(:ubar, p), GA(:ν), Spinor(:v, pp))
    tr_e = fermion_spin_sum(e_M, e_Mc)

    @testset "Individual traces" begin
        # Each trace should have 3 terms with |coeff|=4 and 2 factors
        @test length(tr_mu.terms) == 3
        @test length(tr_e.terms) == 3
        for t in tr_mu.terms
            @test abs(t.coeff) == 4
            @test length(t.factors) == 2
        end
    end

    # ── Step 3: Multiply traces and contract μ,ν ──────────────────
    product = tr_mu * tr_e
    @test length(product.terms) == 9  # 3 × 3

    contracted = contract(product)

    @testset "Contracted result is scalar" begin
        # All Lorentz indices should be contracted away
        for t in contracted.terms
            for f in t.factors
                if f isa Pair
                    @test !(f.a isa LorentzIndex)
                    @test !(f.b isa LorentzIndex)
                end
            end
        end
    end

    # ── Step 4: Substitute Mandelstam variables ───────────────────
    # CM frame kinematics for massless particles:
    #   p·p' = s/2,  k·k' = s/2
    #   p·k = -t/2,  p'·k' = -t/2
    #   p·k' = -u/2, p'·k = -u/2
    #   p² = p'² = k² = k'² = 0
    #
    # Use s=100, θ=π/3: t=-25, u=-75 (s+t+u=0 ✓)

    s_val = 100
    t_val = -25
    u_val = -75

    ctx = SPContext()
    ctx = set_sp(ctx, :p,  :p,  0)          # massless
    ctx = set_sp(ctx, :pp, :pp, 0)
    ctx = set_sp(ctx, :k,  :k,  0)
    ctx = set_sp(ctx, :kp, :kp, 0)
    ctx = set_sp(ctx, :p,  :pp, s_val // 2)       # p·p' = s/2
    ctx = set_sp(ctx, :k,  :kp, s_val // 2)       # k·k' = s/2
    ctx = set_sp(ctx, :p,  :k,  -t_val // 2)      # p·k = -t/2
    ctx = set_sp(ctx, :pp, :kp, -t_val // 2)      # p'·k' = -t/2
    ctx = set_sp(ctx, :p,  :kp, -u_val // 2)      # p·k' = -u/2
    ctx = set_sp(ctx, :pp, :k,  -u_val // 2)      # p'·k = -u/2

    evaluated = evaluate_sp(contracted; ctx=ctx)

    @testset "Numerical evaluation" begin
        @test is_scalar(evaluated)
        # Sum all scalar coefficients
        total = sum(t.coeff for t in evaluated.terms)

        # Expected: 8(t² + u²) = 8(625 + 5625) = 50000
        expected = 8 * (t_val^2 + u_val^2)
        @test total == expected
    end

    # ── Step 5: Verify P&S Eq (5.10) ─────────────────────────────
    @testset "Peskin & Schroeder Eq (5.10)" begin
        # The trace contraction L_μν L'^μν should equal 8(t²+u²)
        # for ANY s,t,u with s+t+u=0.
        # Test with a second set of kinematics: θ=π/4
        s2 = 200
        t2 = -s2 // 2 * (1 - 1 // 2)     # ~ -50 (using cos(π/4)≈1/√2, approximate)
        u2 = -s2 - t2                      # s+t+u=0

        ctx2 = SPContext()
        ctx2 = set_sp(ctx2, :p,  :p,  0)
        ctx2 = set_sp(ctx2, :pp, :pp, 0)
        ctx2 = set_sp(ctx2, :k,  :k,  0)
        ctx2 = set_sp(ctx2, :kp, :kp, 0)
        ctx2 = set_sp(ctx2, :p,  :pp, s2 // 2)
        ctx2 = set_sp(ctx2, :k,  :kp, s2 // 2)
        ctx2 = set_sp(ctx2, :p,  :k,  -t2 // 2)
        ctx2 = set_sp(ctx2, :pp, :kp, -t2 // 2)
        ctx2 = set_sp(ctx2, :p,  :kp, -u2 // 2)
        ctx2 = set_sp(ctx2, :pp, :k,  -u2 // 2)

        eval2 = evaluate_sp(contracted; ctx=ctx2)
        total2 = sum(t.coeff for t in eval2.terms)
        expected2 = 8 * (t2^2 + u2^2)
        @test total2 == expected2

        # The differential cross section:
        # dσ/dΩ = (1/4) Σ|M|² / (64π²s)
        #       = (1/4) × e⁴/s² × L × 1/(64π²s)
        # where L = 8(t²+u²)
        #
        # In CM frame: t²+u² = s²/2 (1+cos²θ)
        # → dσ/dΩ = e⁴ × 2(t²+u²)/(s² × 64π²s)
        #         = e⁴ × 2 × s²/2 × (1+cos²θ) / (s² × 64π²s)
        #         = e⁴(1+cos²θ)/(64π²s)
        #         = (4πα)²(1+cos²θ)/(64π²s)
        #         = α²(1+cos²θ)/(4s)
        #
        # This is P&S Eq (5.10). ✓
        # The algebraic identity is verified by the numerical checks above.
    end
end
