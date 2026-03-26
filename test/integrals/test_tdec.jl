# Tests for Tdec and PaVeReduce
import Feynfeld: Pair

@testset "Tdec rank 0" begin
    result = tdec([], [:p])
    @test result == [(1, 1)]
end

@testset "Tdec rank 1" begin
    result = tdec([(:q, :μ)], [:p])
    @test length(result) == 1
    @test result[1][2] == FVD(:p, :μ)
end

@testset "Tdec rank 2" begin
    result = tdec([(:q, :μ), (:q, :ν)], [:p])
    @test length(result) == 2  # g^{μν} term + p^μ p^ν term
    @test result[1][1] == :C00  # metric coefficient
end

@testset "Tdec rank 2 two external" begin
    result = tdec([(:q, :μ), (:q, :ν)], [:p1, :p2])
    # g^{μν} + p1^μ p1^ν + p1^μ p2^ν + p2^μ p2^ν = 4 terms
    @test length(result) == 4
end

@testset "PaVeReduce: B0 unchanged" begin
    b0 = B0(:pp, :m0, :m1)
    @test pave_reduce(b0) === b0
end

@testset "PaVeReduce: B1 reduces" begin
    b1 = B1(:pp, :m0, :m1)
    result = pave_reduce(b1)
    @test result isa Tuple
    @test result[1] == :B1_reduced
end

@testset "PaVeReduce: B00 reduces" begin
    b00 = B00(:pp, :m0, :m1)
    result = pave_reduce(b00)
    @test result[1] == :B00_reduced
end

@testset "PaVeReduce: A0 unchanged" begin
    a0 = A0(:m)
    @test pave_reduce(a0) === a0
end
