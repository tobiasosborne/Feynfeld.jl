# MUnit translations: FeynCalc ExpandScalarProduct.test
# Source: refs/FeynCalc/Tests/Lorentz/ExpandScalarProduct.test
# Convention: all identities verified symbolically (exact AlgSum equality).
#
# Portable tests: ID1 (FV sum), ID2 (SP sum), ID5 (nested sums with coefficients).
# Skipped: ID3 (selective Momentum filter), ID4 (nested FV — needs coeff*MomentumSum),
# ID6 (mixed D/D-4 dimensions), ID7-9 (LC), ID10-12 (custom tensor/SP assignment),
# ID13-18 (CartesianPair), ID19-21 (selective expansion).

using Test
@isdefined(FeynfeldX) || include(joinpath(@__DIR__, "..", "..", "..", "src", "v2", "FeynfeldX.jl"))
using .FeynfeldX

@testset "MUnit ExpandScalarProduct" begin

    # ---- Helpers ----
    li(s) = LorentzIndex(s, DimD())
    mom(s) = Momentum(s)

    # ==== FourVector expansion ====

    # fcstExpandScalarProduct-ID1
    # Source: refs/FeynCalc/Tests/Lorentz/ExpandScalarProduct.test, ID1
    # Input:  ExpandScalarProduct[FourVector[a + b, mu]]
    # Output: Pair[LorentzIndex[mu], Momentum[a]] + Pair[LorentzIndex[mu], Momentum[b]]
    @testset "ID1: FV[a+b, mu] → FV[a,mu] + FV[b,mu]" begin
        # pair(li(:mu), a+b) creates Pair{LorentzIndex, MomentumSum}
        fv_sum = alg(pair(li(:mu), mom(:a) + mom(:b)))
        result = expand_scalar_product(fv_sum)
        expected = alg(pair(li(:mu), mom(:a))) + alg(pair(li(:mu), mom(:b)))
        @test result == expected
    end

    # ==== Scalar product expansion ====

    # fcstExpandScalarProduct-ID2 (partial: single SP term)
    # Source: refs/FeynCalc/Tests/Lorentz/ExpandScalarProduct.test, ID2
    # "ExpandScalarProduct[SP[a+b, c+d]]" gives 4 terms
    @testset "SP[a+b, c+d] bilinear expansion" begin
        sp_sum = alg(pair(mom(:a) + mom(:b), mom(:c) + mom(:d)))
        result = expand_scalar_product(sp_sum)
        expected = alg(pair(mom(:a), mom(:c))) +
                   alg(pair(mom(:a), mom(:d))) +
                   alg(pair(mom(:b), mom(:c))) +
                   alg(pair(mom(:b), mom(:d)))
        @test result == expected
    end

    # Full ID2: sum of two SP terms with like-term collection
    # SP[a+b, c+d] + SP[a+l, c+d] = 2*SP[a,c] + 2*SP[a,d] + SP[b,c] + SP[b,d] + SP[c,l] + SP[d,l]
    @testset "ID2: SP[a+b,c+d] + SP[a+l,c+d] with collection" begin
        term1 = alg(pair(mom(:a) + mom(:b), mom(:c) + mom(:d)))
        term2 = alg(pair(mom(:a) + mom(:l), mom(:c) + mom(:d)))
        result = expand_scalar_product(term1 + term2)
        expected = 2//1 * alg(pair(mom(:a), mom(:c))) +
                   2//1 * alg(pair(mom(:a), mom(:d))) +
                   alg(pair(mom(:b), mom(:c))) +
                   alg(pair(mom(:b), mom(:d))) +
                   alg(pair(mom(:c), mom(:l))) +
                   alg(pair(mom(:d), mom(:l)))
        @test result == expected
    end

    # ==== Expansion with coefficients ====

    # SP[a-b, c+d] = SP[a,c] + SP[a,d] - SP[b,c] - SP[b,d]
    @testset "SP[a-b, c+d] with minus sign" begin
        sp_diff = alg(pair(mom(:a) - mom(:b), mom(:c) + mom(:d)))
        result = expand_scalar_product(sp_diff)
        expected = alg(pair(mom(:a), mom(:c))) +
                   alg(pair(mom(:a), mom(:d))) -
                   alg(pair(mom(:b), mom(:c))) -
                   alg(pair(mom(:b), mom(:d)))
        @test result == expected
    end

    # SP[2a, b] = 2*SP[a,b]
    @testset "SP[2a, b] = 2*SP[a,b]" begin
        sp_2a = alg(pair(2 * mom(:a), mom(:b)))
        result = expand_scalar_product(sp_2a)
        @test result == 2//1 * alg(pair(mom(:a), mom(:b)))
    end

    # SP[a+b, a+b] = SP[a,a] + 2*SP[a,b] + SP[b,b]
    @testset "SP[a+b, a+b] = a^2 + 2ab + b^2" begin
        sp_sq = alg(pair(mom(:a) + mom(:b), mom(:a) + mom(:b)))
        result = expand_scalar_product(sp_sq)
        expected = alg(pair(mom(:a), mom(:a))) +
                   2//1 * alg(pair(mom(:a), mom(:b))) +
                   alg(pair(mom(:b), mom(:b)))
        @test result == expected
    end

    # SP[p1-p2, p1-p2] = p1^2 - 2*p1.p2 + p2^2
    @testset "SP[p1-p2, p1-p2] Mandelstam-like" begin
        sp_t = alg(pair(mom(:p1) - mom(:p2), mom(:p1) - mom(:p2)))
        result = expand_scalar_product(sp_t)
        expected = alg(pair(mom(:p1), mom(:p1))) -
                   2//1 * alg(pair(mom(:p1), mom(:p2))) +
                   alg(pair(mom(:p2), mom(:p2)))
        @test result == expected
    end

    # ==== FourVector with subtraction ====

    # FV[p1-p2, mu] = FV[p1,mu] - FV[p2,mu]
    @testset "FV[p1-p2, mu] expansion" begin
        fv_diff = alg(pair(li(:mu), mom(:p1) - mom(:p2)))
        result = expand_scalar_product(fv_diff)
        expected = alg(pair(li(:mu), mom(:p1))) - alg(pair(li(:mu), mom(:p2)))
        @test result == expected
    end

    # ==== Product of SP with other factors (algebraic context) ====

    # 3 * SP[a+b, c] = 3*SP[a,c] + 3*SP[b,c]
    @testset "coefficient * SP[a+b, c]" begin
        expr = 3//1 * alg(pair(mom(:a) + mom(:b), mom(:c)))
        result = expand_scalar_product(expr)
        expected = 3//1 * alg(pair(mom(:a), mom(:c))) +
                   3//1 * alg(pair(mom(:b), mom(:c)))
        @test result == expected
    end

    # SP[a+b, c] * MT[mu,nu] — mixed tensor product
    @testset "SP[a+b,c] * MT[mu,nu] product" begin
        expr = alg(pair(mom(:a) + mom(:b), mom(:c))) * alg(pair(li(:mu), li(:nu)))
        result = expand_scalar_product(expr)
        expected = alg(pair(mom(:a), mom(:c))) * alg(pair(li(:mu), li(:nu))) +
                   alg(pair(mom(:b), mom(:c))) * alg(pair(li(:mu), li(:nu)))
        @test result == expected
    end

end
