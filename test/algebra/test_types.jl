# Basic type construction and identity tests for algebra layer
# These validate that the type system is well-formed.

@testset "LorentzIndex" begin
    μ = LorentzIndex(:μ)
    ν = LorentzIndex(:ν)
    @test μ.name == :μ
    @test μ.dim === Dim4()   # default: 4-dimensional (FeynCalc convention)
    @test μ != ν

    μ4 = LorentzIndex(:μ, 4)
    @test μ4.dim === Dim4()

    μD = LorentzIndex(:μ, :D)
    @test μD.dim === DimD()

    # Ordering: alphabetical by name, then by dim
    @test isless(LorentzIndex(:a), LorentzIndex(:b))
    @test isless(LorentzIndex(:μ, Dim4()), LorentzIndex(:μ, DimD()))
end

@testset "FourMomentum (legacy)" begin
    p = FourMomentum(:p)
    k = FourMomentum(:k)
    @test p.name == :p
    @test p != k
    @test p == FourMomentum(:p)
end

@testset "DiracGamma" begin
    μ = LorentzIndex(:μ)
    γμ = DiracGamma(μ)
    @test γμ.index == μ

    γ5 = DiracGamma5()
    @test γ5 isa FeynExpr
end

@testset "Slash notation" begin
    p = FourMomentum(:p)
    pslash = Slash(p)
    @test pslash.momentum == p
end

@testset "Spinors" begin
    p = FourMomentum(:p)
    u = SpinorU(p, :mₑ)
    v = SpinorV(p, :mₑ)
    @test u.momentum == p
    @test v.mass == :mₑ
end

@testset "DiracChain" begin
    p = FourMomentum(:p)
    k = FourMomentum(:k)
    μ = LorentzIndex(:μ)
    chain = DiracChain([SpinorU(p, :m), DiracGamma(μ), SpinorU(k, :m)])
    @test length(chain.elements) == 3
end

@testset "SU(N) colour" begin
    T = SUNMatrix(:T, :a)
    @test T.a == :a

    f = SUNF(:a, :b, :c)
    @test f.a == :a && f.b == :b && f.c == :c

    d = SUND(:a, :b, :c)
    @test d isa FeynExpr

    δ = ColourDelta(:a, :b)
    @test δ.a == :a
end

@testset "PaVe integrals" begin
    a0 = A0(:m)
    @test a0 isa PaVe
    @test length(a0.masses) == 1

    b0 = B0(:p, :m1, :m2)
    @test length(b0.masses) == 2
    @test length(b0.momenta) == 1
end

@testset "Amplitude (Phase 0 placeholder)" begin
    μ = LorentzIndex(:μ)
    t = FTerm(1.0, FeynExpr[DiracGamma(μ)])
    amp = Amplitude([t])
    @test length(amp.terms) == 1
    @test amp.terms[1].coeff == 1.0
end
