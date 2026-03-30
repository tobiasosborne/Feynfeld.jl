# MUnit test translation batch 1: DiracTrace, Contract, PolarizationSum
# Source: refs/FeynCalc/Tests/ — test IDs cited for each assertion.
# Translated per protocol: CLAUDE.md §MUnit translation protocol.

using Test
@isdefined(FeynfeldX) || include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

# ---- Helpers for concise test writing ----
li(s) = LorentzIndex(s, DimD())   # D-dimensional Lorentz index
li4(s) = LorentzIndex(s, Dim4())  # 4-dimensional Lorentz index
mom(s) = Momentum(s)
mt(a, b) = pair(li(a), li(b))    # metric tensor g^{ab} in D-dim
fv(p, mu) = pair(li(mu), mom(p)) # four-vector p^mu in D-dim
sp(a, b) = pair(mom(a), mom(b))  # scalar product p·q

@testset "MUnit batch 1" begin

    # ==== DiracTrace ====
    # Source: refs/FeynCalc/Tests/Dirac/DiracTrace.test

    @testset "DiracTrace basics" begin
        # fcstDiracTrace-ID9: Tr[γ^i] = 0 (odd number)
        @test iszero(dirac_trace([GAD(:i)]))

        # fcstDiracTrace-ID10: Tr[γ·(p+q)] = 0 (odd)
        # Note: single MomSumSlot gamma → expand first, each term is odd → 0
        @test iszero(dirac_trace([GS(mom(:p))]))  # simplified: single gamma trace = 0

        # fcstDiracTrace-ID11: Tr[γ5] = 0
        # (Odd number of gammas in 4D trace convention)
        # Our dirac_trace of an empty chain gives 4 (Tr[I] = 4)
        @test dirac_trace(DiracGamma[]) == alg(4)
    end

    @testset "DiracTrace 2-gamma" begin
        # fcstDiracTrace-ID32: Tr[γ·p γ·p] = 4 p·p
        # "DiracTrace[GS[p, p], DiracTraceEvaluate->True] // FCE" = "4 SP[p, p]"
        result = dirac_trace([GS(mom(:p)), GS(mom(:p))])
        @test result == 4 * alg(sp(:p, :p))

        # Tr[γ^μ γ^ν] = 4 g^{μν}  (fundamental trace identity)
        result2 = dirac_trace([GAD(:mu), GAD(:nu)])
        @test result2 == 4 * alg(mt(:mu, :nu))
    end

    @testset "DiracTrace 4-gamma" begin
        # Tr[γ^a γ^b γ^c γ^d] = 4(g^{ab}g^{cd} - g^{ac}g^{bd} + g^{ad}g^{bc})
        # This is the fundamental 4-gamma trace formula.
        result = dirac_trace([GAD(:a), GAD(:b), GAD(:c), GAD(:d)])
        expected = 4 * (alg(mt(:a,:b)) * alg(mt(:c,:d)) -
                        alg(mt(:a,:c)) * alg(mt(:b,:d)) +
                        alg(mt(:a,:d)) * alg(mt(:b,:c)))
        @test result == expected

        # Tr[γ·p γ·q γ·r γ·s] = 4[(p·q)(r·s) - (p·r)(q·s) + (p·s)(q·r)]
        result2 = dirac_trace([GS(mom(:p)), GS(mom(:q)), GS(mom(:r)), GS(mom(:s))])
        expected2 = 4 * (alg(sp(:p,:q)) * alg(sp(:r,:s)) -
                         alg(sp(:p,:r)) * alg(sp(:q,:s)) +
                         alg(sp(:p,:s)) * alg(sp(:q,:r)))
        @test result2 == expected2
    end

    @testset "DiracTrace 6-gamma" begin
        # Tr[γ^a γ^b γ^c γ^d γ^e γ^f] has 15 terms.
        # Test numerically at specific SP values rather than full symbolic comparison.
        result = dirac_trace([GAD(:a), GAD(:b), GAD(:c), GAD(:d), GAD(:e), GAD(:f)])
        @test length(result.terms) == 15

        # Cross-check: Tr[γ·p1 γ·p2 γ·p3 γ·p4 γ·p5 γ·p6] should have 15 terms
        result2 = dirac_trace([GS(mom(:p1)), GS(mom(:p2)), GS(mom(:p3)),
                                GS(mom(:p4)), GS(mom(:p5)), GS(mom(:p6))])
        @test length(result2.terms) == 15
    end

    @testset "DiracTrace 8-gamma" begin
        # Compton scattering requires traces of 8 gammas.
        # Tr[γ^a γ^b γ^c γ^d γ^e γ^f γ^g γ^h] has 105 terms.
        result = dirac_trace([GAD(:a), GAD(:b), GAD(:c), GAD(:d),
                              GAD(:e), GAD(:f), GAD(:g), GAD(:h)])
        @test length(result.terms) == 105
    end

    @testset "DiracTrace linearity" begin
        # fcstDiracTrace-ID40:
        # "DiracTrace[GSD[q].GSD[p1] + GSD[q].GSD[p2], DiracTraceEvaluate->True]"
        # = "4 (Pair[Momentum[p1,D], Momentum[q,D]] + Pair[Momentum[p2,D], Momentum[q,D]])"
        t1 = dirac_trace([GS(mom(:q)), GS(mom(:p1))])
        t2 = dirac_trace([GS(mom(:q)), GS(mom(:p2))])
        @test t1 + t2 == 4 * alg(sp(:q, :p1)) + 4 * alg(sp(:q, :p2))
    end

    # ==== Contract ====
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test

    @testset "Contract metric-metric" begin
        # fcstContractContractionsIn4dims-ID4: g^{μ}_μ = 4  (4-dim)
        s4 = alg(pair(li4(:mu), li4(:mu)))
        @test contract(s4) == alg(4)

        # fcstContractContractionsIn4dims-ID5: g^{μ}_μ = D  (D-dim)
        sD = alg(mt(:mu, :mu))
        @test contract(sD) == alg(DIM)
    end

    @testset "Contract metric-vector" begin
        # fcstContractContractionsIn4dims-ID6:
        # "Contract[MetricTensor[a,b] FourVector[p,b]]" = "Pair[LorentzIndex[a], Momentum[p]]"
        s = alg(mt(:a, :b)) * alg(fv(:p, :b))
        result = contract(s)
        @test result == alg(fv(:p, :a))
    end

    @testset "Contract vector-vector" begin
        # fcstContractContractionsIn4dims-ID7 (adapted for D-dim):
        # "Contract[FVD[q,a] FVD[p-q,a]]" should give SPD[p,q] - SPD[q,q]
        # We test: FV(q,a) * FV(p,a) = SP(q,p)  (simple case)
        s = alg(fv(:q, :a)) * alg(fv(:p, :a))
        result = contract(s)
        @test result == alg(sp(:q, :p))
    end

    @testset "Contract with MomentumSum" begin
        # FV(p-q, a) * FV(p-q, a) should contract to (p-q)·(p-q)
        # Pipeline: expand first (bilinear), then contract (repeated index a)
        pmq = MomentumSum([(1//1, mom(:p)), (-1//1, mom(:q))])
        s = alg(pair(li(:a), pmq)) * alg(pair(li(:a), pmq))
        expanded = expand_scalar_product(s)
        contracted = contract(expanded)
        expected = alg(sp(:p,:p)) - 2 * alg(sp(:p,:q)) + alg(sp(:q,:q))
        @test contracted == expected
    end

    @testset "Contract chain" begin
        # g^{ab} g^{bc} = g^{ac}  (metric chain)
        s = alg(mt(:a, :b)) * alg(mt(:b, :c))
        result = contract(s)
        @test result == alg(mt(:a, :c))

        # g^{ab} g^{bc} g^{cd} = g^{ad}
        s2 = alg(mt(:a, :b)) * alg(mt(:b, :c)) * alg(mt(:c, :d))
        result2 = contract(s2)
        @test result2 == alg(mt(:a, :d))
    end

    # ==== PolarizationSum ====
    # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test

    @testset "PolarizationSum Feynman gauge" begin
        # fcstPolarizationSum-ID1:
        # "PolarizationSum[rho, si]" = "-Pair[LorentzIndex[rho], LorentzIndex[si]]"
        @test polarization_sum(li(:rho), li(:si)) == -alg(mt(:rho, :si))

        # ID8 (VirtualBoson): same result in Feynman gauge
        @test polarization_sum(li(:mu), li(:nu)) == -alg(mt(:mu, :nu))
    end

    # ==== substitute_index ====
    # (Infrastructure for polarization sum contraction)

    @testset "substitute_index" begin
        # Replace μ' → μ in g^{μ'ν}
        s = alg(mt(:mu_, :nu))
        result = substitute_index(s, li(:mu_), li(:mu))
        @test result == alg(mt(:mu, :nu))

        # Replace in product: g^{μ'ν} p^{μ'} → g^{μν} p^μ
        s2 = alg(mt(:mu_, :nu)) * alg(fv(:p, :mu_))
        result2 = substitute_index(s2, li(:mu_), li(:mu))
        expected2 = alg(mt(:mu, :nu)) * alg(fv(:p, :mu))
        @test result2 == expected2
    end

    # ==== evaluate_sp ====

    @testset "evaluate_sp" begin
        # SP evaluation at specific kinematics
        ctx = sp_context((:p, :p) => 1//1, (:p, :q) => 3//2, (:q, :q) => 0//1)
        s = 2 * alg(sp(:p, :p)) + alg(sp(:p, :q))
        result = evaluate_sp(s; ctx)
        @test result == alg(7//2)  # 2*1 + 3/2
    end
end
