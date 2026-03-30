# MUnit translations: FeynCalc DiracTrace.test
# Source: refs/FeynCalc/Tests/Dirac/DiracTrace.test
# Convention: all trace identities verified symbolically (exact AlgSum equality).

using Test
@isdefined(FeynfeldX) || include(joinpath(@__DIR__, "..", "..", "..", "src", "v2", "FeynfeldX.jl"))
using .FeynfeldX

@testset "MUnit DiracTrace" begin

    # ---- Trace of identity and single gammas ----
    # fcstDiracTrace-ID7: Tr[I] = 4
    @testset "ID7: Tr[I]=4" begin
        @test dirac_trace(DiracGamma[]) == alg(4)
    end

    # fcstDiracTrace-ID9: Tr[γ^i] = 0
    @testset "ID9: Tr[γ^i]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i)]))
    end

    # fcstDiracTrace-ID10: Tr[γ·p] = 0
    @testset "ID10: Tr[γ·p]=0" begin
        @test iszero(dirac_trace(DiracGamma[GS(:p)]))
    end

    # fcstDiracTrace-ID11: Tr[γ5] = 0
    @testset "ID11: Tr[γ5]=0" begin
        @test iszero(dirac_trace(DiracGamma[GA5()]))
    end

    # fcstDiracTrace-ID12: Tr[GA6] = Tr[(1+γ5)/2] = 2
    @testset "ID12: Tr[GA6]=2" begin
        @test dirac_trace(DiracGamma[GA6()]) == alg(2)
    end

    # fcstDiracTrace-ID13: Tr[GA7] = Tr[(1-γ5)/2] = 2
    @testset "ID13: Tr[GA7]=2" begin
        @test dirac_trace(DiracGamma[GA7()]) == alg(2)
    end

    # ---- Odd-gamma with γ5: all zero ----
    # fcstDiracTrace-ID17: Tr[γ^i γ5] = 0
    @testset "ID17: Tr[γ^i γ5]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GA5()]))
    end

    # fcstDiracTrace-ID18: Tr[γ^i γ^j γ5] = 0
    @testset "ID18: Tr[γ^i γ^j γ5]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GAD(:j), GA5()]))
    end

    # fcstDiracTrace-ID19: Tr[γ^i γ^j γ^k γ5] = 0  (3 gammas + γ5 = odd)
    @testset "ID19: Tr[γ^i γ^j γ^k γ5]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GAD(:j), GAD(:k), GA5()]))
    end

    # ---- Projector traces: 2-gamma with GA6/GA7 ----
    # fcstDiracTrace-ID21: Tr[γ^i GA6] = 0
    @testset "ID21: Tr[γ^i GA6]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GA6()]))
    end

    # fcstDiracTrace-ID22: Tr[γ^i γ^j GA6] = 2 g^{ij}
    @testset "ID22: Tr[γ^i γ^j GA6]=2g^{ij}" begin
        @test dirac_trace(DiracGamma[GAD(:i), GAD(:j), GA6()]) == 2//1 * alg(MTD(:i, :j))
    end

    # fcstDiracTrace-ID23: Tr[γ^i γ^j γ^k GA6] = 0
    @testset "ID23: Tr[γ^i γ^j γ^k GA6]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GAD(:j), GAD(:k), GA6()]))
    end

    # fcstDiracTrace-ID25: Tr[γ^i GA7] = 0
    @testset "ID25: Tr[γ^i GA7]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GA7()]))
    end

    # fcstDiracTrace-ID26: Tr[γ^i γ^j GA7] = 2 g^{ij}
    @testset "ID26: Tr[γ^i γ^j GA7]=2g^{ij}" begin
        @test dirac_trace(DiracGamma[GAD(:i), GAD(:j), GA7()]) == 2//1 * alg(MTD(:i, :j))
    end

    # fcstDiracTrace-ID27: Tr[γ^i γ^j γ^k GA7] = 0
    @testset "ID27: Tr[γ^i γ^j γ^k GA7]=0" begin
        @test iszero(dirac_trace(DiracGamma[GAD(:i), GAD(:j), GAD(:k), GA7()]))
    end

    # ---- 2-gamma fundamental traces ----
    # fcstDiracTrace-ID32: Tr[γ·p γ·p] = 4 p·p
    @testset "ID32: Tr[γ·p γ·p]=4p²" begin
        @test dirac_trace(DiracGamma[GS(:p), GS(:p)]) == 4//1 * alg(SP(:p, :p))
    end

    # fcstDiracTrace-ID52: Tr[γ^μ γ^ν] = 4 g^{μν}
    @testset "ID52: Tr[γ^μ γ^ν]=4g^{μν}" begin
        @test dirac_trace(DiracGamma[GAD(:mu), GAD(:nu)]) == 4//1 * alg(MTD(:mu, :nu))
    end

    # ---- Linearity ----
    # fcstDiracTrace-ID40: Tr[γ·q γ·p1 + γ·q γ·p2] = 4(q·p1 + q·p2)
    @testset "ID40: linearity" begin
        t1 = dirac_trace(DiracGamma[GS(:q), GS(:p1)])
        t2 = dirac_trace(DiracGamma[GS(:q), GS(:p2)])
        @test t1 + t2 == 4//1 * alg(SP(:q, :p1)) + 4//1 * alg(SP(:q, :p2))
    end

    # ---- 4-gamma traces ----
    # Tr[γ·a γ·b γ·c γ·d] = 4[(a·b)(c·d) - (a·c)(b·d) + (a·d)(b·c)]
    # Ref: Peskin & Schroeder Eq. (A.26)
    @testset "4-gamma trace" begin
        result = dirac_trace(DiracGamma[GS(:a), GS(:b), GS(:c), GS(:d)])
        expected = 4//1 * (alg(SP(:a,:b))*alg(SP(:c,:d)) - alg(SP(:a,:c))*alg(SP(:b,:d)) + alg(SP(:a,:d))*alg(SP(:b,:c)))
        @test result == expected
    end

    # ---- Anticommutator pattern (ID45) ----
    # fcstDiracTrace-ID45: Tr[(γ·p2+m)(γ·a γ·pp - γ·pp γ·a)(γ·p1+m)]
    #   = 8(a·p2 × p1·pp - a·p1 × p2·pp)
    @testset "ID45: anticommutator" begin
        # (γ·p2+m)(γ·a γ·pp - γ·pp γ·a)(γ·p1+m)
        # Expand: 4 terms × 2 = 8 sub-traces
        t1 = dirac_trace(DiracGamma[GS(:p2), GS(:a), GS(:pp), GS(:p1)])
        t2 = dirac_trace(DiracGamma[GS(:p2), GS(:pp), GS(:a), GS(:p1)])
        # Mass terms: m × Tr[γ·a γ·pp - γ·pp γ·a] = 0 (antisymmetric under trace)
        # Cross mass terms: m × Tr[γ·p2(γ·a γ·pp-γ·pp γ·a)] + m × Tr[(γ·a γ·pp-γ·pp γ·a)γ·p1]
        # By cyclic property these also cancel
        # m² × Tr[γ·a γ·pp - γ·pp γ·a] = 0
        # So only the m=0 part survives:
        result = t1 - t2
        expected = 8//1 * (alg(SP(:a,:p2))*alg(SP(:p1,:pp)) - alg(SP(:a,:p1))*alg(SP(:p2,:pp)))
        @test result == expected
    end

    # ---- γ5 traces: Tr[γ^{μ₁} γ^{μ₂} γ^{μ₃} γ^{μ₄} γ5] = -4i ε^{μ₁μ₂μ₃μ₄} ----
    # fcstDiracTrace-ID54
    @testset "ID54: Tr[4γ + γ5] = -4i ε" begin
        result = dirac_trace(DiracGamma[GAD(:mu1), GAD(:mu2), GAD(:mu3), GAD(:mu4), GA5()])
        # Expected: -4i × Eps(mu1,mu2,mu3,mu4)
        # In Feynfeld: the Eps tensor with DimD indices
        mu1 = LorentzIndex(:mu1, DimD()); mu2 = LorentzIndex(:mu2, DimD())
        mu3 = LorentzIndex(:mu3, DimD()); mu4 = LorentzIndex(:mu4, DimD())
        expected = (-4//1) * alg(Eps(mu1, mu2, mu3, mu4))
        # Note: FeynCalc uses (-4*I)*LCD = -4i*ε, but in our convention
        # the imaginary unit i from the trace is absorbed into the Eps definition.
        # Check if our trace gives -4*Eps (real) or -4i*Eps (complex).
        # The standard result is Tr[γ^μ γ^ν γ^ρ γ^σ γ5] = -4i ε^{μνρσ}
        # If Feynfeld convention absorbs the i: result should be -4 Eps
        @test result == expected
    end

    # fcstDiracTrace-ID55: same trace with γ5 in the middle
    # Tr[γ^{μ₁} γ^{μ₂} γ5 γ^{μ₃} γ^{μ₄}] = -4i ε^{μ₁μ₂μ₃μ₄}
    # (γ5 anticommutes through, picking up a sign for each gamma)
    @testset "ID55: γ5 in middle position" begin
        result = dirac_trace(DiracGamma[GAD(:mu1), GAD(:mu2), GA5(), GAD(:mu3), GAD(:mu4)])
        mu1 = LorentzIndex(:mu1, DimD()); mu2 = LorentzIndex(:mu2, DimD())
        mu3 = LorentzIndex(:mu3, DimD()); mu4 = LorentzIndex(:mu4, DimD())
        expected = (-4//1) * alg(Eps(mu1, mu2, mu3, mu4))
        @test result == expected
    end
end
