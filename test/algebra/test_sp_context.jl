# Tests for SPContext (scalar product value storage)
# Ref: FeynCalc ScalarProduct.m, FCClearScalarProducts.m

@testset "Basic set/get" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :p, :m²)
    @test get_sp(ctx, :p, :p) == :m²
    @test get_sp(ctx, :k, :k) === nothing
end

@testset "Symmetry" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :q, :s)
    @test get_sp(ctx, :p, :q) == :s
    @test get_sp(ctx, :q, :p) == :s  # symmetric
end

@testset "Immutability" begin
    ctx1 = SPContext()
    ctx2 = set_sp(ctx1, :p, :p, :m²)
    # Original unchanged
    @test get_sp(ctx1, :p, :p) === nothing
    @test get_sp(ctx2, :p, :p) == :m²
end

@testset "Multiple assignments" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :p, :m²)
    ctx = set_sp(ctx, :p, :q, :s)
    ctx = set_sp(ctx, :q, :q, 0)
    @test get_sp(ctx, :p, :p) == :m²
    @test get_sp(ctx, :p, :q) == :s
    @test get_sp(ctx, :q, :q) == 0
end

@testset "Override" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :q, :old)
    ctx = set_sp(ctx, :p, :q, :new)
    @test get_sp(ctx, :p, :q) == :new
end
