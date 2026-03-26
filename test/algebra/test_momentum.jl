# Tests for Momentum type
# Ref: FeynCalc SharedObjects.m Momentum[p, dim]

@testset "Construction" begin
    p = Momentum(:p)
    @test p.name == :p
    @test p.dim === Dim4()  # default: 4D

    pD = Momentum(:p, DimD())
    @test pD.dim === DimD()

    pE = Momentum(:p, DimDm4())
    @test pE.dim === DimDm4()

    # Convenience: integer/symbol dim
    @test Momentum(:p, 4).dim === Dim4()
    @test Momentum(:p, :D).dim === DimD()
end

@testset "Equality" begin
    @test Momentum(:p) == Momentum(:p)
    @test Momentum(:p) != Momentum(:q)
    @test Momentum(:p, Dim4()) != Momentum(:p, DimD())
end

@testset "Ordering" begin
    # Alphabetical by name
    @test isless(Momentum(:a), Momentum(:b))
    @test isless(Momentum(:k), Momentum(:p))
    @test !isless(Momentum(:p), Momentum(:k))

    # Same name: order by dim (Dim4 < DimD < DimDm4)
    @test isless(Momentum(:p, Dim4()), Momentum(:p, DimD()))
    @test isless(Momentum(:p, DimD()), Momentum(:p, DimDm4()))
end
