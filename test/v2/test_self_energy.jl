# 1-loop electron self-energy in QED — cockroach detector
#
# Σ(p) numerator: γ^μ (m + γ·q) γ_μ  where q = p-k
# After D-dim contraction: D·m·I + (2-D)·γ·q
#
# Tests DiracExpr (matrix-valued), DiracTrick (gamma contraction),
# MomentumSum expansion in traces, and DimPoly coefficients.

using Test
@isdefined(FeynfeldX) || include(joinpath(@__DIR__, "..", "..", "src", "v2", "FeynfeldX.jl"))
using .FeynfeldX

@testset "1-loop self-energy" begin

    @testset "DiracExpr construction" begin
        # Can build matrix-valued expressions
        de_scalar = DiracExpr(alg(3))
        @test !isempty(de_scalar.terms)

        de_gamma = DiracExpr(DiracChain([GS(:p)]))
        @test length(de_gamma.terms) == 1

        # Arithmetic
        sum = de_scalar + de_gamma
        @test length(sum.terms) == 2

        # Scalar * DiracExpr
        scaled = DIM * de_gamma
        @test length(scaled.terms) == 1
    end

    @testset "DiracTrick n=0: γ^μ γ_μ = D" begin
        mu = LorentzIndex(:mu, DimD())
        chain = DiracChain([DiracGamma(LISlot(mu)), DiracGamma(LISlot(mu))])
        de = DiracExpr(chain)
        result = dirac_trick(de)

        # Should be D * identity
        @test length(result.terms) == 1
        coeff, ch = result.terms[1]
        @test isempty(ch.elements)  # identity matrix
        @test haskey(coeff.terms, FactorKey())
        @test coeff.terms[FactorKey()] == DIM
    end

    @testset "DiracTrick n=1: γ^μ γ·p γ_μ = -(D-2) γ·p" begin
        mu = LorentzIndex(:mu, DimD())
        p = Momentum(:p, DimD())
        chain = DiracChain([DiracGamma(LISlot(mu)), GS(p), DiracGamma(LISlot(mu))])
        de = DiracExpr(chain)
        result = dirac_trick(de)

        # Should be (2-D) * γ·p
        @test length(result.terms) == 1
        coeff, ch = result.terms[1]
        @test length(ch.elements) == 1  # single gamma
        @test haskey(coeff.terms, FactorKey())
        val = coeff.terms[FactorKey()]
        # 2 - D as DimPoly
        @test val == (2 - DIM) || val == DimPoly([2, -1])
    end

    @testset "Self-energy numerator: γ^μ (m + γ·q) γ_μ" begin
        mu = LorentzIndex(:mu, DimD())
        q = Momentum(:q, DimD())  # q = p - k
        m = 1//1  # mass (using rational for exactness)

        # Build: γ^μ (m·I + γ·q) γ_μ as DiracExpr
        g_mu = DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))
        scalar_part = DiracExpr(m * alg(1))  # m * identity
        vector_part = DiracExpr(DiracChain([GS(q)]))  # γ·q
        inner = scalar_part + vector_part

        # Full numerator: γ^μ · (m + γ·q) · γ_μ
        numerator = g_mu * inner * g_mu
        result = dirac_trick(numerator)

        # Expected: D*m*I + (2-D)*γ·q
        # So we should have exactly 2 terms
        @test length(result.terms) == 2

        # Find the scalar (identity) and vector (γ·q) terms
        scalar_term = nothing
        vector_term = nothing
        for (coeff, ch) in result.terms
            if isempty(ch.elements)
                scalar_term = coeff
            elseif length(ch.elements) == 1
                vector_term = coeff
            end
        end

        @test scalar_term !== nothing
        @test vector_term !== nothing

        if scalar_term !== nothing
            # Scalar coefficient should be D*m = D (since m=1)
            @test haskey(scalar_term.terms, FactorKey())
            @test scalar_term.terms[FactorKey()] == DIM
        end

        if vector_term !== nothing
            # Vector coefficient should be (2-D)
            @test haskey(vector_term.terms, FactorKey())
            @test vector_term.terms[FactorKey()] == (2 - DIM)
        end

        # Evaluate at D=4: should give 4m*I - 2*γ·q
        if scalar_term !== nothing && vector_term !== nothing
            s_val = evaluate_dim(scalar_term.terms[FactorKey()])
            v_val = evaluate_dim(vector_term.terms[FactorKey()])
            @test s_val == 4   # 4m (m=1)
            @test v_val == -2  # -(4-2) = -2
        end
    end

    @testset "MomentumSum in trace: Tr[γ·(p-k) γ·q]" begin
        p = Momentum(:p)
        k = Momentum(:k)
        pk = p - k  # MomentumSum
        @test pk isa MomentumSum

        g_pk = DiracGamma(MomSumSlot(pk))
        tr = dirac_trace([g_pk, GS(:q)])

        # Expected: Tr[γ·p γ·q] - Tr[γ·k γ·q] = 4(p·q) - 4(k·q)
        @test !iszero(tr)

        # Should have two terms: 4*SP(p,q) and -4*SP(k,q)
        @test length(tr.terms) == 2
    end

    @testset "Self-energy trace check: Tr[Σ(p) γ·p]" begin
        # Taking trace of (D*m*I + (2-D)*γ·q) with an extra γ·p:
        # Tr[D*m*γ·p] + Tr[(2-D)*γ·q*γ·p]
        # = D*m*0 (odd trace) + (2-D)*4*(q·p)
        # = 4(2-D)(q·p)
        #
        # At D=4: -8(q·p)

        q = Momentum(:q, DimD())
        p = Momentum(:p, DimD())

        # Tr[(2-D) γ·q γ·p] — the only surviving term
        tr = dirac_trace([GS(q), GS(p)])
        # = 4 q·p (in D-dim, Tr[γ·a γ·b] = 4 a·b since trace of 2 gammas = 4*metric)
        @test !iszero(tr)
    end
end
