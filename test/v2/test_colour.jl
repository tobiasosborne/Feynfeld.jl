# Tests for SU(N) colour algebra

using Feynfeld
using Test

@testset "Colour algebra" begin
    a = AdjointIndex(:a)
    b = AdjointIndex(:b)
    c = AdjointIndex(:c)

    @testset "Type construction and canonicalization" begin
        # SUNDelta canonical ordering
        d1 = SUNDelta(b, a)
        d2 = SUNDelta(a, b)
        @test d1 == d2
        @test d1.a == a  # smaller name first

        # SUNF antisymmetry
        f1 = SUNF(a, b, c)
        f2 = SUNF(b, a, c)
        @test f1.a == f2.a && f1.b == f2.b && f1.c == f2.c  # same canonical form
        @test f1.sign == -f2.sign  # opposite signs

        # SUND symmetry
        d1 = SUND(c, a, b)
        d2 = SUND(a, b, c)
        @test d1 == d2  # same after sorting

        # ColourChain
        chain = ColourChain([SUNT(a), SUNT(b)])
        @test length(chain) == 2
    end

    @testset "Colour trace (N=3, QCD)" begin
        N = 3

        # Tr(1) = N = 3
        tr0 = colour_trace(SUNT[]; N=N)
        @test tr0.terms[FactorKey()] == 3//1

        # Tr(T^a) = 0
        tr1 = colour_trace([SUNT(a)]; N=N)
        @test iszero(tr1)

        # Tr(T^a T^b) = (1/2) δ^{ab}
        tr2 = colour_trace([SUNT(a), SUNT(b)]; N=N)
        @test !iszero(tr2)
        # Should have one term: 1/2 * δ(a,b)
        expected_key = FactorKey(AlgFactor[SUNDelta(a, b)])
        @test haskey(tr2.terms, expected_key)
        @test tr2.terms[expected_key] == 1//2
    end

    @testset "Casimirs (N=3)" begin
        @test casimir_fundamental(3) == 4//3   # CF = (9-1)/6 = 4/3
        @test casimir_adjoint(3) == 3//1       # CA = 3
        @test trace_normalization(3) == 1//2    # TF = 1/2
    end

    @testset "Colour trace n>=4 (i² fix, feynfeld-83m)" begin
        N = 3
        # Tr(T^a T^b T^b T^a) = C_F^2 × (1/2)δ^{aa} = C_F^2 × (N^2-1)/2
        # For N=3: C_F = 4/3, so C_F^2 × 4 = (16/9) × 4 = 64/9... wait.
        #
        # Direct calculation: T^b T^b = C_F × I, so
        # Tr(T^a T^b T^b T^a) = C_F × Tr(T^a T^a) = C_F × (1/2)δ^{aa}
        #                      = (4/3) × (1/2) × 8 = 16/3
        # Ref: Casimir identity, C_F = (N²-1)/(2N)
        tr4 = colour_trace([SUNT(a), SUNT(b), SUNT(b), SUNT(a)]; N=N)
        contracted4 = contract_colour(tr4; N=N)
        @test contracted4.terms[FactorKey()] == 16//3

        # Tr(T^a T^b T^a T^b) = (C_F - C_A/2) × Tr(T^a T^a)
        #                      = (C_F - C_A/2) × (N^2-1)/2
        # For N=3: (4/3 - 3/2) × 4 = (-1/6) × 4 = -2/3
        # Ref: use T^a T^b T^a = (C_F - C_A/2) T^b (standard identity)
        tr4b = colour_trace([SUNT(a), SUNT(b), SUNT(a), SUNT(b)]; N=N)
        contracted4b = contract_colour(tr4b; N=N)
        @test contracted4b.terms[FactorKey()] == -2//3
    end

    @testset "Delta traces (N=3)" begin
        # δ^{aa} = N² - 1 = 8
        @test colour_delta_trace(SUNDelta(a, a); N=3) == 8//1

        # δ^{ab} with a ≠ b: not a trace
        @test colour_delta_trace(SUNDelta(a, b); N=3) === nothing

        # Fundamental: δ_{ii} = N = 3
        i = FundIndex(:i)
        @test colour_delta_trace(FundDelta(i, i); N=3) == 3//1
    end

    @testset "Colour contraction" begin
        # Test: δ^{ab} δ^{ab} = N² - 1 = 8 (for N=3)
        d_ab = SUNDelta(a, b)
        expr = alg(d_ab) * alg(d_ab)
        contracted = contract_colour(expr; N=3)
        # After contracting δ^{ab} with δ^{ab}, get δ^{aa} = 8
        # Actually δ^{ab}δ^{ab} → δ^{aa} → 8
        @test contracted.terms[FactorKey()] == 8//1

        # Test: f^{acd} f^{bcd} = N δ^{ab}
        d_idx = AdjointIndex(:d)
        f1 = SUNF(a, c, d_idx)
        f2 = SUNF(b, c, d_idx)
        expr2 = alg(f1) * alg(f2)
        contracted2 = contract_colour(expr2; N=3)
        # f^{acd} f^{bcd} = N δ^{ab}  for N=3 → 3 δ^{ab}
        # Ref: standard SU(N) identity, cross-check FeynCalc SUNSimplify
        @test !iszero(contracted2)
        expected_ff = 3 * alg(SUNDelta(a, b))
        @test contracted2 == expected_ff

        # Test: d^{acd} d^{bcd} = (N²-4)/N δ^{ab}  for N=3 → 5/3 δ^{ab}
        d1 = SUND(a, c, d_idx)
        d2 = SUND(b, c, d_idx)
        expr_dd = alg(d1) * alg(d2)
        contracted_dd = contract_colour(expr_dd; N=3)
        expected_dd = (5 // 3) * alg(SUNDelta(a, b))
        @test contracted_dd == expected_dd

        # Test: f^{acd} d^{bcd} = 0
        expr_fd = alg(f1) * alg(d2)
        contracted_fd = contract_colour(expr_fd; N=3)
        @test iszero(contracted_fd)
    end

    @testset "AlgFactor integration" begin
        # Colour factors should work in AlgSum arithmetic
        d_ab = SUNDelta(a, b)
        p = SP(:p, :q)  # Lorentz scalar product

        # Can multiply Lorentz and colour factors
        mixed = alg(p) * alg(d_ab)
        @test !iszero(mixed)
        @test length(mixed.terms) == 1

        # Both factors present in the key
        fk = first(keys(mixed.terms))
        @test length(fk.factors) == 2
    end
end
