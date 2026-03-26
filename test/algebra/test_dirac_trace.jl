# Tests for DiracTrace
# Ref: FeynCalc Tests/Dirac/DiracTrace.test
import Feynfeld: Pair

@testset "Trace base cases" begin
    # Tr[1] = 4
    @test dirac_trace(DiracChain(DiracElement[])) == 4

    # Tr[single gamma] = 0 (odd)
    @test dirac_trace(dot(GA(:μ))) == 0

    # Tr[3 gammas] = 0 (odd)
    @test dirac_trace(dot(GA(:μ), GA(:ν), GA(:ρ))) == 0

    # Tr[gamma5] = 0
    @test dirac_trace(dot(GA5())) == 0
end

@testset "Trace of 2 gammas" begin
    # Tr[g^μ g^ν] = 4 g^{μν}
    # Ref: DiracTrace.test standard trace formula
    result = dirac_trace(dot(GA(:μ), GA(:ν)))
    expected = Feynfeld._mul_coeff(4, MT(:μ, :ν))
    @test result == expected
end

@testset "Trace of 4 gammas" begin
    # Tr[g^μ g^ν g^ρ g^σ] = 4(g^{μν}g^{ρσ} - g^{μρ}g^{νσ} + g^{μσ}g^{νρ})
    result = dirac_trace(dot(GA(:μ), GA(:ν), GA(:ρ), GA(:σ)))
    # Result is a list of 3 terms
    @test result isa Vector
    @test length(result) == 3
end

@testset "Trace with slashed momenta" begin
    # Tr[p-slash q-slash] = 4 p·q
    result = dirac_trace(dot(GS(:p), GS(:q)))
    expected = Feynfeld._mul_coeff(4, SP(:p, :q))
    @test result == expected

    # Tr[p-slash] = 0 (odd)
    @test dirac_trace(dot(GS(:p))) == 0
end

@testset "Trace with gamma5 and 4 gammas" begin
    # Tr[g5 g^a g^b g^c g^d] = -4 * ε^{abcd} (with LCS = -1)
    result = dirac_trace(dot(GA5(), GA(:a), GA(:b), GA(:c), GA(:d)))
    @test result != 0
    # Result should involve Eps
    @test result isa Expr || result isa Eps
end

@testset "Trace gamma5 base cases" begin
    # Tr[g5 g^a g^b] = 0 (only 2 regular gammas, need 4)
    @test dirac_trace(dot(GA5(), GA(:a), GA(:b))) == 0

    # Tr[g5 g5] = Tr[1] = 4 (even g5 count)
    @test dirac_trace(dot(GA5(), GA5())) == 4

    # Tr[g5 g5 g^a g^b] = Tr[g^a g^b] = 4 g^{ab}
    result = dirac_trace(dot(GA5(), GA5(), GA(:a), GA(:b)))
    @test result == Feynfeld._mul_coeff(4, MT(:a, :b))
end

@testset "Trace spinor rejection" begin
    u = Spinor(:u, Momentum(:p), :m)
    @test_throws ErrorException dirac_trace(dot(u, GA(:μ)))
end
