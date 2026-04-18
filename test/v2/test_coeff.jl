# Tests for DimPoly coefficient algebra

using Feynfeld
using Test

@testset "DimPoly" begin
    @testset "Construction and display" begin
        @test DimPoly([4]) == DimPoly([4//1])
        @test DimPoly([0, 1]) == DIM
        @test DimPoly([-4, 1]) == DIM_MINUS_4
        @test iszero(DimPoly(Int[]))
        @test isone(DimPoly([1]))
        @test string(DIM) == "D"
        s = string(DIM_MINUS_4)
        @test s == "D - 4" || s == "-4 + D"  # accept either display order
    end

    @testset "Arithmetic" begin
        @test DIM + DIM == DimPoly([0, 2])              # 2D
        @test DIM - DIM_MINUS_4 == DimPoly([4])          # D - (D-4) = 4
        @test 3 * DIM == DimPoly([0, 3])                 # 3D
        @test DIM * DIM == DimPoly([0, 0, 1])            # D^2 ← SymDim would error here!
        @test DIM_MINUS_4 * DIM == DimPoly([0, -4, 1])   # D(D-4) = D^2 - 4D
        @test (DIM - 4) * (DIM + 4) == DimPoly([-16, 0, 1])  # (D-4)(D+4) = D^2-16
    end

    @testset "Evaluate at D=4" begin
        @test evaluate_dim(DIM) == 4
        @test evaluate_dim(DIM_MINUS_4) == 0
        @test evaluate_dim(DIM * DIM) == 16        # D^2 at D=4
        @test evaluate_dim(DimPoly([-16, 0, 1])) == 0  # D^2-16 at D=4
        @test evaluate_dim(42) == 42
    end

    @testset "Normalisation" begin
        @test normalise_coeff(DimPoly([7//1])) === 7//1   # constant DimPoly → Rational
        @test normalise_coeff(DimPoly(Int[])) === 0//1     # zero DimPoly → 0
        @test normalise_coeff(DIM) === DIM                 # non-constant stays DimPoly
        @test normalise_coeff(5) === 5//1
    end

    @testset "mul_coeff / add_coeff" begin
        @test mul_coeff(3//1, 4//1) === 12//1
        @test mul_coeff(DIM, 2//1) == DimPoly([0, 2])
        @test mul_coeff(DIM, DIM) == DimPoly([0, 0, 1])  # D*D = D^2
        @test add_coeff(DIM, 4//1) == DimPoly([4, 1])    # D+4
        @test add_coeff(DIM, DIM) == DimPoly([0, 2])      # D+D = 2D
    end

    @testset "Integration with AlgSum" begin
        # DimPoly coefficients should flow through AlgSum arithmetic
        s1 = alg(DIM)  # AlgSum with scalar D
        s2 = alg(4)    # AlgSum with scalar 4
        sum = s1 + s2
        @test !iszero(sum)

        # D * AlgSum should work
        p = SP(:p, :q)
        sp_sum = alg(p)
        scaled = DIM * sp_sum
        @test !iszero(scaled)
    end
end
