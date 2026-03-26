# Tests for SU(N) colour algebra
# Ref: FeynCalc Tests/SUN/SUNTrace.test, SUNSimplify.test

@testset "SUNIndex ordering" begin
    @test isless(SUNIndex(:a), SUNIndex(:b))
    @test SUNDelta(SUNIndex(:b), SUNIndex(:a)).a == SUNIndex(:a)
end

@testset "SUNF antisymmetry" begin
    a, b, c = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c)
    f1 = SUNF(a, b, c)
    f2 = SUNF(b, a, c)
    # Same indices but different sign
    @test f1.a == f2.a && f1.b == f2.b && f1.c == f2.c
    @test f1.sign == -f2.sign

    # Repeated index → structure constants would be zero in contraction
    # (SUNF itself doesn't enforce this; SUNSimplify does)
end

@testset "SUND symmetry" begin
    a, b, c = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c)
    d1 = SUND(a, b, c)
    d2 = SUND(c, a, b)
    @test d1 == d2  # symmetric, same canonical form
end

@testset "ColourChain" begin
    a, b = SUNIndex(:a), SUNIndex(:b)
    chain = ColourChain([SUNT(a), SUNT(b)])
    @test length(chain) == 2
end

# ── SUNTrace ─────────────────────────────────────────────────────────

@testset "SUNTrace: Tr(1) = N" begin
    @test sun_trace(ColourChain(SUNT[])) == :N
end

@testset "SUNTrace: Tr(T^a) = 0" begin
    @test sun_trace(SUNT(SUNIndex(:a))) == 0
end

@testset "SUNTrace: Tr(T^a T^b) = (1/2) δ^{ab}" begin
    a, b = SUNIndex(:a), SUNIndex(:b)
    result = sun_trace(SUNT(a), SUNT(b))
    @test result isa Tuple
    @test result[1] == 1 // 2
    @test result[2] == SUNDelta(a, b)
end

@testset "SUNTrace: Tr(T^a T^b T^c) = (1/4)(d + if)" begin
    a, b, c = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c)
    result = sun_trace(SUNT(a), SUNT(b), SUNT(c))
    @test result isa Tuple
    @test result[1] == 1 // 4
    @test result[2] isa SUND
    @test result[3] isa SUNF
end

@testset "SUNTrace: Tr(T^a T^b T^c T^d) recursive" begin
    a, b, c, d = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c), SUNIndex(:d)
    result = sun_trace(SUNT(a), SUNT(b), SUNT(c), SUNT(d))
    # Returns a structured tuple with recursive sub-traces
    @test result isa Tuple
    @test length(result) == 5  # (delta, sub1, d, f, sub2)
end

# ── Delta trace ──────────────────────────────────────────────────────

@testset "Delta trace" begin
    a = SUNIndex(:a)
    i = SUNFIndex(:i)
    # δ^{aa} = N²-1
    @test delta_trace(SUNDelta(a, a)) == :(N ^ 2 - 1)
    # δ_F^{ii} = N
    @test delta_trace(SUNFDelta(i, i)) == :N
    # Non-diagonal unchanged
    @test delta_trace(SUNDelta(SUNIndex(:a), SUNIndex(:b))) isa SUNDelta
end

# ── Structure constant contractions ──────────────────────────────────

@testset "f^{acd} f^{bcd} = N δ^{ab}" begin
    a, b, c, d = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c), SUNIndex(:d)
    f1 = SUNF(a, c, d)
    f2 = SUNF(b, c, d)
    result = contract_ff(f1, f2)
    @test result !== nothing
    @test result[1] == :N
    @test result[2] == SUNDelta(a, b)
end

@testset "f^{abc} d^{abd} = 0" begin
    a, b, c, d = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c), SUNIndex(:d)
    f = SUNF(a, b, c)
    dd = SUND(a, b, d)
    @test contract_fd(f, dd) == 0
end

@testset "f^{abc} f^{abc} = 2 CA² CF" begin
    a, b, c = SUNIndex(:a), SUNIndex(:b), SUNIndex(:c)
    f1 = SUNF(a, b, c)
    f2 = SUNF(a, b, c)
    result = contract_ff_full(f1, f2)
    @test result !== nothing
    @test result isa Expr
end
