# Tests for PaVe{N} types and numerical evaluation.
# Ground truth: textbook formulas + known special values.

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "PaVe type system" begin
    @testset "Construction" begin
        @test A0(1.0) isa PaVe{1}
        @test B0(4.0, 0.0, 1.0) isa PaVe{2}
        @test B1(4.0, 0.0, 1.0) isa PaVe{2}
        @test C0(1.0, 2.0, 3.0, 0.1, 0.2, 0.3) isa PaVe{3}

        @test A0(1.0).masses == [1.0]
        @test A0(1.0).invariants == Float64[]
        @test B0(4.0, 0.0, 1.0).invariants == [4.0]
        @test B0(4.0, 0.0, 1.0).masses == [0.0, 1.0]
        @test B1(4.0, 0.0, 1.0).indices == [1]
    end

    @testset "Validation" begin
        # Wrong number of invariants
        @test_throws ArgumentError PaVe{2}(Int[], [1.0, 2.0], [0.0, 1.0])
        # Wrong number of masses
        @test_throws ArgumentError PaVe{2}(Int[], [1.0], [0.0])
        # N=0 not allowed
        @test_throws ArgumentError PaVe{0}(Int[], Float64[], Float64[])
    end

    @testset "Equality and hashing" begin
        @test B0(4.0, 0.0, 1.0) == B0(4.0, 0.0, 1.0)
        @test B0(4.0, 0.0, 1.0) != B0(4.0, 0.0, 2.0)
        @test A0(1.0) != B0(1.0, 0.0, 0.0)  # different N
        @test hash(B0(4.0, 0.0, 1.0)) == hash(B0(4.0, 0.0, 1.0))
        @test hash(A0(1.0)) != hash(A0(2.0))

        # Dict key usage
        d = Dict(B0(4.0, 0.0, 1.0) => 1, A0(1.0) => 2)
        @test d[B0(4.0, 0.0, 1.0)] == 1
        @test d[A0(1.0)] == 2
    end

    @testset "Index sorting" begin
        pv = PaVe{2}([2, 1], [1.0], [0.0, 1.0])
        @test pv.indices == [1, 2]  # sorted
    end

    @testset "Show" begin
        @test repr(A0(1.0)) == "A(1.0)"
        @test repr(B0(4.0, 0.0, 1.0)) == "B(4.0, 0.0, 1.0)"
        @test repr(B1(4.0, 0.0, 1.0)) == "B1(4.0, 0.0, 1.0)"
    end
end

@testset "A0 numerical evaluation" begin
    # A0(0) = 0 (massless tadpole vanishes in dim-reg)
    @test evaluate(A0(0.0)) == 0.0 + 0.0im

    # A0(m²) = m²(1 - ln(m²/μ²)) at μ² = 1
    # A0(1.0) = 1.0 * (1 - ln(1)) = 1.0
    @test real(evaluate(A0(1.0))) ≈ 1.0 atol = 1e-14

    # A0(e) at μ² = 1: e * (1 - ln(e)) = e * (1 - 1) = 0
    @test real(evaluate(A0(exp(1.0)))) ≈ 0.0 atol = 1e-12

    # A0(m²) at μ² = m²: m²(1 - ln(1)) = m²
    @test real(evaluate(A0(5.0); mu2 = 5.0)) ≈ 5.0 atol = 1e-14
end

@testset "B0 numerical evaluation" begin
    @testset "Both masses zero" begin
        # B0(p², 0, 0) = 2 - ln(-p²/μ²) for p² < 0 (spacelike)
        # At p² = -1.0, μ² = 1: B0 = 2 - ln(1) = 2
        @test real(evaluate(B0(-1.0, 0.0, 0.0))) ≈ 2.0 atol = 1e-10
        @test imag(evaluate(B0(-1.0, 0.0, 0.0))) ≈ 0.0 atol = 1e-10

        # At p² = -e, μ² = 1: B0 = 2 - ln(e) = 1
        @test real(evaluate(B0(-exp(1.0), 0.0, 0.0))) ≈ 1.0 atol = 1e-10

        # For p² > 0 (timelike): imaginary part = π
        @test imag(evaluate(B0(1.0, 0.0, 0.0))) ≈ π atol = 1e-10
        @test real(evaluate(B0(1.0, 0.0, 0.0))) ≈ 2.0 atol = 1e-10
    end

    @testset "One mass zero" begin
        # B0(0, 0, m²) = 1 - ln(m²/μ²)
        # At m² = 1, μ² = 1: B0 = 1 - 0 = 1
        @test real(evaluate(B0(0.0, 0.0, 1.0))) ≈ 1.0 atol = 1e-10

        # At m² = e, μ² = 1: B0 = 1 - 1 = 0
        @test real(evaluate(B0(0.0, 0.0, exp(1.0)))) ≈ 0.0 atol = 1e-10

        # B0(m², 0, m²) = 2 - ln(m²/μ²) (p² = m² special case)
        @test real(evaluate(B0(1.0, 0.0, 1.0))) ≈ 2.0 atol = 1e-10
    end

    @testset "Equal masses" begin
        # B0(0, m², m²) = -ln(m²/μ²)
        # At m² = 1, μ² = 1: B0 = 0
        @test real(evaluate(B0(0.0, 1.0, 1.0))) ≈ 0.0 atol = 1e-10

        # At m² = e, μ² = 1: B0 = -1
        @test real(evaluate(B0(0.0, exp(1.0), exp(1.0)))) ≈ -1.0 atol = 1e-10

        # B0(4m², m², m²) = 2 - ln(m²/μ²) (threshold, β=0)
        # Near-threshold has mild integrand peak; quadrature tolerance ~1e-7
        m2 = 1.0
        @test real(evaluate(B0(4 * m2, m2, m2))) ≈ 2.0 - log(m2) atol = 1e-6

        # B0 with nonzero p² below threshold (p² < 4m²): should be real
        @test abs(imag(evaluate(B0(1.0, 1.0, 1.0)))) < 1e-10

        # B0 above threshold (p² > 4m²): Im(B₀) = πβ where β = √(1-4m²/p²)
        # Ref: derived from Denner 1993 Eq. (4.23), -iε prescription
        # For p²=10, m²=1: β = √(1-4/10) = √(3/5), Im = π√(3/5)
        @test imag(evaluate(B0(10.0, 1.0, 1.0))) ≈ π * sqrt(3 / 5) atol = 1e-10

        # B0 above threshold, unequal masses: Im = π√λ/p²
        # λ(p², m₀², m₁²) = p⁴ - 2p²(m₀²+m₁²) + (m₀²-m₁²)²
        # For p²=20, m₀²=1, m₁²=4: λ = 400 - 2(20)(5) + 9 = 209
        lambda_uneq = 20.0^2 - 2 * 20 * (1.0 + 4.0) + (1.0 - 4.0)^2
        @test imag(evaluate(B0(20.0, 1.0, 4.0))) ≈ π * sqrt(lambda_uneq) / 20.0 atol = 1e-8

        # One massless above threshold: Im = π(p²-m²)/p²
        # For p²=5, m₀²=0, m₁²=1: Im = π(5-1)/5 = 4π/5
        @test imag(evaluate(B0(5.0, 0.0, 1.0))) ≈ 4π / 5 atol = 1e-10

        # Verify B0(p², m², m²) = 2 - ln(m²/μ²) - β ln((β+1)/(β-1))
        # at p² = -4 (spacelike), m² = 1: β² = 1 + 1 = 2, β = √2
        p2, m2_test = -4.0, 1.0
        β = sqrt(1.0 - 4.0 * m2_test / p2)  # real for spacelike
        expected_mm = 2.0 - log(m2_test) - β * log((β + 1) / (β - 1))
        @test real(evaluate(B0(p2, m2_test, m2_test))) ≈ expected_mm atol = 1e-10
    end

    @testset "Zero momentum, different masses" begin
        # B0(0, m₁², m₂²) = 1 - ln(m₂²/μ²) + m₁²/(m₁²-m₂²) ln(m₂²/m₁²)
        m02, m12 = 1.0, 4.0
        expected = 1.0 - log(m12) + m02 / (m02 - m12) * log(m12 / m02)
        @test real(evaluate(B0(0.0, m02, m12))) ≈ expected atol = 1e-10
    end

    @testset "B0 mass symmetry" begin
        # B0(p², m₁², m₂²) = B0(p², m₂², m₁²) by Feynman parameter x→(1-x)
        p2 = 5.0
        b0_a = evaluate(B0(p2, 0.0, 2.0))
        b0_b = evaluate(B0(p2, 2.0, 0.0))
        @test real(b0_a) ≈ real(b0_b) atol = 1e-10
        @test imag(b0_a) ≈ imag(b0_b) atol = 1e-10

        # Also test with two different nonzero masses
        b0_c = evaluate(B0(p2, 1.0, 3.0))
        b0_d = evaluate(B0(p2, 3.0, 1.0))
        @test real(b0_c) ≈ real(b0_d) atol = 1e-10
        @test imag(b0_c) ≈ imag(b0_d) atol = 1e-10
    end

    @testset "General B0 (two different nonzero masses)" begin
        # B0 at spacelike momentum: should be real
        b0 = evaluate(B0(-4.0, 1.0, 2.0))
        @test abs(imag(b0)) < 1e-10

        # Consistency: general formula at equal masses should match equal-mass formula
        b0_gen = evaluate(B0(5.0, 1.0, 1.0 + 1e-14))  # slightly perturbed
        b0_eq = evaluate(B0(5.0, 1.0, 1.0))            # exact equal
        @test real(b0_gen) ≈ real(b0_eq) atol = 1e-6
        @test imag(b0_gen) ≈ imag(b0_eq) atol = 1e-6
    end
end

@testset "B1 numerical evaluation" begin
    # B1 = [A0(m₁²) - A0(m₀²) - (p²+m₀²-m₁²)B0] / (2p²)
    # Test at specific kinematics
    p2, m02, m12 = 4.0, 1.0, 1.0
    b1 = evaluate(B1(p2, m02, m12))
    # Verify via manual computation
    a0_m1 = evaluate(A0(m12))
    a0_m0 = evaluate(A0(m02))
    b0 = evaluate(B0(p2, m02, m12))
    expected = (a0_m1 - a0_m0 - (p2 + m02 - m12) * b0) / (2 * p2)
    @test real(b1) ≈ real(expected) atol = 1e-10
    @test imag(b1) ≈ imag(expected) atol = 1e-10
end
