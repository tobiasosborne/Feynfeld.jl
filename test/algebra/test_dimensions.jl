# Tests for BMHV dimensional regularisation algebra
# Validates the complete BMHV projection rules from FeynCalc.
# Ref: FeynCalc SharedObjects.m, PairContract.m dimEval rules

@testset "Dimension types" begin
    @test Dim4() isa Dim4
    @test DimD() isa DimD
    @test DimDm4() isa DimDm4
end

@testset "BMHV projection algebra" begin
    # Same-dim projections (identity)
    @test dim_contract(Dim4(), Dim4()) isa Dim4
    @test dim_contract(DimD(), DimD()) isa DimD
    @test dim_contract(DimDm4(), DimDm4()) isa DimDm4

    # 4D ∩ D = 4D (D projected onto 4 gives 4)
    @test dim_contract(Dim4(), DimD()) isa Dim4
    @test dim_contract(DimD(), Dim4()) isa Dim4

    # D ∩ (D-4) = D-4 (D projected onto D-4 gives D-4)
    @test dim_contract(DimD(), DimDm4()) isa DimDm4
    @test dim_contract(DimDm4(), DimD()) isa DimDm4

    # 4 ∩ (D-4) = 0 (vanishes)
    @test dim_contract(Dim4(), DimDm4()) === nothing
    @test dim_contract(DimDm4(), Dim4()) === nothing
end

@testset "Dimension traces" begin
    # g^μ_μ = 4 in 4D
    @test dim_trace(Dim4()) == 4
    # g^μ_μ = D in D dimensions
    @test dim_trace(DimD()) == :D
    # g^μ_μ = D-4 in evanescent
    @test dim_trace(DimDm4()) == :(D - 4)
end

@testset "to_dim conversion" begin
    @test to_dim(4) isa Dim4
    @test to_dim(:D) isa DimD
    @test to_dim(:(D - 4)) isa DimDm4

    # Idempotent
    @test to_dim(Dim4()) isa Dim4
    @test to_dim(DimD()) isa DimD
    @test to_dim(DimDm4()) isa DimDm4

    # Invalid inputs
    @test_throws ErrorException to_dim(3)
    @test_throws ErrorException to_dim(:E)
end
