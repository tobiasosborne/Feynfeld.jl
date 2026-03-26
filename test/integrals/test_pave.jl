# Tests for PaVe type system and FeynAmpDenominator
# Ref: FeynCalc Tests/LoopIntegrals/PaVe.test

# ── PaVe construction ────────────────────────────────────────────────

@testset "PaVe{1}: A-functions" begin
    a0 = A0(:m²)
    @test a0 isa PaVe{1}
    @test a0.indices == Int[]
    @test a0.masses == Any[:m²]
    @test a0.invariants == Any[]

    a00 = A00(:m²)
    @test a00.indices == [0, 0]
end

@testset "PaVe{2}: B-functions" begin
    b0 = B0(:pp, :m0², :m1²)
    @test b0 isa PaVe{2}
    @test b0.indices == Int[]
    @test b0.invariants == Any[:pp]

    b1 = B1(:pp, :m0², :m1²)
    @test b1.indices == [1]

    b00 = B00(:pp, :m0², :m1²)
    @test b00.indices == [0, 0]

    b11 = B11(:pp, :m0², :m1²)
    @test b11.indices == [1, 1]
end

@testset "PaVe{3}: C-functions" begin
    c0 = C0(:p10, :p12, :p20, :m0², :m1², :m2²)
    @test c0 isa PaVe{3}
    @test length(c0.invariants) == 3
    @test length(c0.masses) == 3

    # Tensor C-function
    c12 = PaVe{3}([1, 2], Any[:p10, :p12, :p20], Any[:m0², :m1², :m2²])
    @test c12.indices == [1, 2]  # sorted
end

@testset "PaVe{4}: D-functions" begin
    d0 = D0(:p10, :p12, :p23, :p30, :p20, :p13, :m0², :m1², :m2², :m3²)
    @test d0 isa PaVe{4}
    @test length(d0.invariants) == 6
    @test length(d0.masses) == 4
end

@testset "PaVe index sorting" begin
    # fcstPaVe-ID12: PaVe[2,1,1,...] → PaVe[1,1,2,...]
    pave = PaVe{3}([2, 1, 1], Any[:p10, :p12, :p20], Any[:m0, :m1, :m2])
    @test pave.indices == [1, 1, 2]
end

@testset "PaVe validation" begin
    # Wrong number of invariants
    @test_throws ErrorException PaVe{2}(Int[], Any[:p1, :p2], Any[:m0, :m1])
    # Wrong number of masses
    @test_throws ErrorException PaVe{2}(Int[], Any[:pp], Any[:m0])
end

@testset "PaVe equality" begin
    @test B0(:pp, :m0, :m1) == B0(:pp, :m0, :m1)
    @test A0(:m) != A0(:M)
    @test hash(B0(:pp, :m0, :m1)) == hash(B0(:pp, :m0, :m1))
end

# ── FeynAmpDenominator ──────────────────────────────────────────────

@testset "FAD construction" begin
    # Massless propagator
    fad = FAD(:q)
    @test fad isa FeynAmpDenominator
    @test length(fad) == 1

    # Massive propagator
    fad = FAD((:q, :m))
    @test fad.propagators[1].mass == :m

    # Multi-propagator
    fad = FAD(:q, (:q, :m))
    @test length(fad) == 2
end

@testset "FAD multiplication" begin
    f1 = FAD(:q)
    f2 = FAD((:q, :m))
    f3 = f1 * f2
    @test length(f3) == 2
end

@testset "Propagator types" begin
    pd = PropagatorDenominator(Momentum(:q), :m)
    @test pd.momentum == Momentum(:q)

    sp = StandardPropagator(Momentum(:q), 0, :m², 1, 1)
    @test sp.power == 1

    gp = GenericPropagator(:(q^2 - m^2))
    @test gp.power == 1
end
