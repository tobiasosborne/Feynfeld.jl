# Tracer bullet: e+e- → mu+mu- at tree level
# Validates Peskin & Schroeder Eq (5.10): |M|^2 = 8(t^2 + u^2) (massless, unpolarised)
#
# This is the same physics as test/test_ee_mumu.jl but using FeynfeldX's
# Julia-idiomatic design: parametric types, dispatch, SymDim coefficients.

using Test

# Load FeynfeldX module
@isdefined(FeynfeldX) || include(joinpath(@__DIR__, "..", "..", "src", "v2", "FeynfeldX.jl"))
using .FeynfeldX

@testset "FeynfeldX: e+e- → mu+mu- tracer bullet" begin

    # ---- Step 1: Define momenta ----
    p1 = Momentum(:p1)  # incoming e-
    p2 = Momentum(:p2)  # incoming e+
    k1 = Momentum(:k1)  # outgoing mu-
    k2 = Momentum(:k2)  # outgoing mu+

    @testset "Momentum arithmetic" begin
        @test p1 + p2 isa MomentumSum
        @test p1 - p1 === nothing  # cancels to zero
    end

    # ---- Step 2: Build amplitude chains ----
    # QED tree-level: e+e- → gamma* → mu+mu-
    # M = (-ie)^2 / s * [vbar(p2) γ^μ u(p1)] [ubar(k1) γ_μ v(k2)]
    # For |M|^2 we need two chains and their conjugates.

    @testset "Dirac chain construction" begin
        mu = LorentzIndex(:mu)
        chain1 = dot(vbar(p2), GA(:mu), u(p1))
        chain2 = dot(ubar(k1), GA(:mu), v(k2))
        @test chain1 isa DiracChain
        @test chain2 isa DiracChain
        @test length(chain1.elements) == 3
    end

    # ---- Step 3: Spin sum via completeness relations ----
    @testset "Spin sum (massless)" begin
        # For |M|^2, after spin-summing we get traces.
        # Chain 1: vbar(p2) γ^mu u(p1) → conjugate gives ubar(p1) γ^nu v(p2)
        # |M_1|^2 → Tr[p2-slash γ^mu p1-slash γ^nu]  (massless)
        # |M_2|^2 → Tr[k1-slash γ_mu k2-slash γ_nu]

        # Build the trace chains directly (as the spin sum would produce)
        mu_idx = LorentzIndex(:mu)
        nu_idx = LorentzIndex(:nu)

        # Trace 1: Tr[p2-slash γ^mu p1-slash γ^nu]
        trace1_gammas = [GS(:p2), GA(:mu), GS(:p1), GA(:nu)]
        tr1 = dirac_trace(trace1_gammas)
        @test tr1 isa AlgSum
        @test !iszero(tr1)

        # Trace 2: Tr[k1-slash γ_mu k2-slash γ_nu]
        trace2_gammas = [GS(:k1), GA(:mu), GS(:k2), GA(:nu)]
        tr2 = dirac_trace(trace2_gammas)
        @test tr2 isa AlgSum
        @test !iszero(tr2)
    end

    # ---- Step 4: Compute traces and contract ----
    @testset "Trace computation" begin
        # Tr[a-slash γ^mu b-slash γ^nu] = 4(a^mu b^nu + a^nu b^mu - (a.b) g^{mu nu})
        trace1_gammas = [GS(:p2), GA(:mu), GS(:p1), GA(:nu)]
        tr1 = dirac_trace(trace1_gammas)

        trace2_gammas = [GS(:k1), GA(:mu), GS(:k2), GA(:nu)]
        tr2 = dirac_trace(trace2_gammas)

        # |M|^2 ∝ tr1 * tr2 (contracted over mu, nu)
        product = tr1 * tr2
        @test product isa AlgSum

        # Contract Lorentz indices mu, nu
        contracted = contract(product)
        @test contracted isa AlgSum
    end

    # ---- Step 5: Full pipeline → Mandelstam → P&S (5.10) ----
    @testset "Full pipeline: P&S Eq (5.10)" begin
        # Compute traces
        tr1 = dirac_trace([GS(:p2), GA(:mu), GS(:p1), GA(:nu)])
        tr2 = dirac_trace([GS(:k1), GA(:mu), GS(:k2), GA(:nu)])

        # Multiply and contract
        product = tr1 * tr2
        contracted = contract(product)

        # Expand scalar products (distribute MomentumSum if any)
        expanded = expand_scalar_product(contracted)

        # Assign Mandelstam variables:
        # s = (p1+p2)^2, t = (p1-k1)^2, u = (p1-k2)^2
        # Massless: p1^2 = p2^2 = k1^2 = k2^2 = 0
        # From these: p1.p2 = s/2, k1.k2 = s/2,
        #             p1.k1 = -t/2, p2.k2 = -t/2,
        #             p1.k2 = -u/2, p2.k1 = -u/2
        #
        # For the final check we evaluate numerically at a specific kinematic point.
        # Choose: s=10, t=-3, u=-7 (s+t+u=0 for massless)

        s_val = 10//1
        t_val = -3//1
        u_val = -7//1

        ctx = sp_context(
            (:p1, :p1) => 0//1,
            (:p2, :p2) => 0//1,
            (:k1, :k1) => 0//1,
            (:k2, :k2) => 0//1,
            (:p1, :p2) => s_val // 2,
            (:k1, :k2) => s_val // 2,
            (:p1, :k1) => -t_val // 2,
            (:p2, :k2) => -t_val // 2,
            (:p1, :k2) => -u_val // 2,
            (:p2, :k1) => -u_val // 2,
        )

        result = evaluate_sp(expanded; ctx=ctx)

        # Result should be a pure number (all Lorentz indices contracted,
        # all scalar products evaluated)
        @test length(result.terms) == 1
        scalar_key = FactorKey()
        @test haskey(result.terms, scalar_key)

        numerical = result.terms[scalar_key]

        # P&S Eq (5.10): 8(t^2 + u^2) (with e=1 coupling stripped)
        expected = 8 * (t_val^2 + u_val^2)
        @test numerical == expected

        # Print the result for human verification
        println("  e+e- → mu+mu- |M|^2 = $(numerical)")
        println("  Expected 8(t²+u²)   = $(expected)")
        println("  t=$t_val, u=$u_val → 8($(t_val^2) + $(u_val^2)) = $(expected)")
    end
end
