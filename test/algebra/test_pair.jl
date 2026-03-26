# Tests for Pair type and convenience constructors
# Ref: FeynCalc SharedObjects.m Pair[x, y]
import Feynfeld: Pair

@testset "Pair construction" begin
    # Metric tensor: Pair(LorentzIndex, LorentzIndex)
    g = Pair(LorentzIndex(:μ), LorentzIndex(:ν))
    @test g isa Pair
    @test g.a == LorentzIndex(:μ)
    @test g.b == LorentzIndex(:ν)

    # Scalar product: Pair(Momentum, Momentum)
    sp = Pair(Momentum(:p), Momentum(:q))
    @test sp.a == Momentum(:p)
    @test sp.b == Momentum(:q)

    # Four-vector: Pair(LorentzIndex, Momentum)
    fv = Pair(LorentzIndex(:μ), Momentum(:p))
    @test fv.a isa LorentzIndex
    @test fv.b isa Momentum
end

@testset "Canonical ordering (Orderless)" begin
    # LorentzIndex comes before Momentum
    p1 = Pair(Momentum(:p), LorentzIndex(:μ))
    p2 = Pair(LorentzIndex(:μ), Momentum(:p))
    @test p1 == p2
    @test p1.a isa LorentzIndex  # LorentzIndex is slot a

    # Alphabetical within same type
    sp = Pair(Momentum(:q), Momentum(:p))
    @test sp.a == Momentum(:p)  # p < q
    @test sp.b == Momentum(:q)

    g = Pair(LorentzIndex(:ν), LorentzIndex(:μ))
    @test g.a == LorentzIndex(:μ)  # μ < ν
    @test g.b == LorentzIndex(:ν)
end

@testset "BMHV projection in Pair" begin
    # D meets 4 → projects to 4
    p = Pair(LorentzIndex(:μ, DimD()), LorentzIndex(:ν, Dim4()))
    @test p.a.dim === Dim4()
    @test p.b.dim === Dim4()

    # D meets (D-4) → projects to (D-4)
    p = Pair(Momentum(:p, DimD()), Momentum(:q, DimDm4()))
    @test p.a.dim === DimDm4()
    @test p.b.dim === DimDm4()

    # Same dim → no change
    p = Pair(LorentzIndex(:μ, DimD()), LorentzIndex(:ν, DimD()))
    @test p.a.dim === DimD()
end

@testset "BMHV vanishing: pair() factory" begin
    # 4 meets (D-4) → vanishes
    @test pair(LorentzIndex(:μ, Dim4()), LorentzIndex(:ν, DimDm4())) == 0
    @test pair(Momentum(:p, Dim4()), Momentum(:q, DimDm4())) == 0

    # Non-vanishing cases return Pair
    @test pair(LorentzIndex(:μ, Dim4()), LorentzIndex(:ν, Dim4())) isa Pair
    @test pair(LorentzIndex(:μ, DimD()), LorentzIndex(:ν, Dim4())) isa Pair

    # Direct Pair() errors on vanishing
    @test_throws ErrorException Pair(LorentzIndex(:μ, Dim4()), LorentzIndex(:ν, DimDm4()))
end

@testset "Convenience: MT/MTD/MTE" begin
    g = MT(:μ, :ν)
    @test g isa Pair
    @test g.a.dim === Dim4()

    gD = MTD(:μ, :ν)
    @test gD.a.dim === DimD()

    gE = MTE(:μ, :ν)
    @test gE.a.dim === DimDm4()
end

@testset "Convenience: FV/FVD/FVE" begin
    fv = FV(:p, :μ)
    @test fv isa Pair
    @test fv.a == LorentzIndex(:μ)
    @test fv.b == Momentum(:p)

    fvD = FVD(:p, :μ)
    @test fvD.a.dim === DimD()
    @test fvD.b.dim === DimD()

    fvE = FVE(:p, :μ)
    @test fvE.a.dim === DimDm4()
end

@testset "Convenience: SP/SPD/SPE" begin
    sp = SP(:p, :q)
    @test sp isa Pair
    @test sp.a == Momentum(:p)
    @test sp.b == Momentum(:q)

    # p² shorthand
    sp2 = SP(:p)
    @test sp2.a == Momentum(:p)
    @test sp2.b == Momentum(:p)

    spD = SPD(:p, :q)
    @test spD.a.dim === DimD()

    spE = SPE(:p, :q)
    @test spE.a.dim === DimDm4()
end

@testset "Equality and hashing" begin
    p1 = SP(:p, :q)
    p2 = SP(:q, :p)  # reversed → same canonical form
    @test p1 == p2
    @test hash(p1) == hash(p2)

    # Use as Dict key
    d = Dict(p1 => 42)
    @test d[p2] == 42

    # LorentzIndex and Momentum with same name don't hash-collide
    li = LorentzIndex(:p)
    m = Momentum(:p)
    @test hash(li) != hash(m)
end

@testset "Edge cases" begin
    # Diagonal metric g^{μ,μ} — valid before trace evaluation
    g_diag = Pair(LorentzIndex(:μ), LorentzIndex(:μ))
    @test g_diag.a == g_diag.b

    # Same-dim evanescent pair (non-vanishing)
    @test pair(LorentzIndex(:μ, DimDm4()), LorentzIndex(:ν, DimDm4())) isa Pair

    # Same-dim D-dimensional pair
    @test pair(Momentum(:p, DimD()), Momentum(:q, DimD())) isa Pair
end
