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
    γμ = GA(:μ)
    @test γμ isa DiracGamma
    @test γμ.slot isa LISlot
    @test γμ.slot.index == LorentzIndex(:μ)

    γ5 = GA5()
    @test γ5 isa DiracGamma
    @test γ5.slot isa SpecialSlot
    @test γ5.slot.id == 5

    # Slashed momentum
    pslash = GS(:p)
    @test pslash.slot isa MomSlot
    @test pslash.slot.mom == Momentum(:p)

    # D-dimensional
    γμD = GAD(:μ)
    @test gamma_dim(γμD) === DimD()

    # BMHV vanishing: 4D index in (D-4) slot
    @test dirac_gamma(LorentzIndex(:μ), DimDm4()) == 0
end

@testset "Spinor" begin
    u = Spinor(:u, Momentum(:p), :mₑ)
    v = Spinor(:v, Momentum(:p), :mₑ)
    @test u.kind == :u
    @test u.momentum == Momentum(:p)
    @test v.mass == :mₑ
    @test_throws ErrorException Spinor(:x, Momentum(:p))
end

@testset "DiracChain" begin
    chain = DiracChain([GA(:μ), GA(:ν)])
    @test length(chain) == 2
end

@testset "SU(N) colour types" begin
    a = SUNIndex(:a)
    b = SUNIndex(:b)
    i = SUNFIndex(:i)
    j = SUNFIndex(:j)

    # Generator
    T = SUNT(a)
    @test T.a == a

    # Structure constants
    f = SUNF(a, b, SUNIndex(:c))
    @test f isa SUNF
    @test f.sign == 1  # already sorted

    # Antisymmetry: f(b,a,c) has sign -1
    f2 = SUNF(b, a, SUNIndex(:c))
    @test f2.sign == -1

    # Symmetric d
    d = SUND(a, b, SUNIndex(:c))
    @test d isa SUND

    # Deltas (canonical ordering)
    δ = SUNDelta(b, a)
    @test δ.a == a  # sorted: a < b

    δF = SUNFDelta(j, i)
    @test δF.i == i  # sorted: i < j

    # SUNTF
    tf = SUNTF([a, b], i, j)
    @test length(tf.adj) == 2
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
    t = FTerm(1.0, FeynExpr[GA(:μ)])
    amp = Amplitude([t])
    @test length(amp.terms) == 1
    @test amp.terms[1].coeff == 1.0
end
