# Tests for ExpandScalarProduct
# Ref: FeynCalc Tests/Lorentz/ExpandScalarProduct.test
import Feynfeld: Pair

# ── MomentumSum construction ─────────────────────────────────────────

@testset "MomentumSum arithmetic" begin
    p = Momentum(:p)
    q = Momentum(:q)

    # Addition
    pq = p + q
    @test pq isa MomentumSum
    @test length(pq.terms) == 2

    # Subtraction
    pmq = p - q
    @test pmq isa MomentumSum

    # Self-subtraction → nothing (zero)
    @test (p - p) === nothing

    # Self-addition → MomentumSum with one term, coeff 2
    pp = p + p
    @test pp isa MomentumSum
    @test length(pp.terms) == 1
    @test pp.terms[1] == (2 // 1, :p)

    # Scalar multiplication
    @test (2 * p) isa MomentumSum
    @test (2 * p).terms == [(2 // 1, :p)]

    # Negation
    neg_p = -p
    @test neg_p isa MomentumSum
    @test neg_p.terms == [(-1 // 1, :p)]

    # Three-term sum
    k = Momentum(:k)
    pqk = p + q + k
    @test pqk isa MomentumSum
    @test length(pqk.terms) == 3
end

@testset "MomentumSum dimension check" begin
    pD = Momentum(:p, DimD())
    q4 = Momentum(:q, Dim4())
    @test_throws ErrorException pD + q4
end

# ── Pair with MomentumSum ────────────────────────────────────────────

@testset "Pair with MomentumSum" begin
    pq = Momentum(:p) + Momentum(:q)

    # FourVector with sum: (p+q)^μ
    fv = Pair(LorentzIndex(:μ), pq)
    @test fv isa Pair
    @test fv.b isa MomentumSum

    # Scalar product with sum: (p+q)·k
    sp = Pair(Momentum(:k), pq)
    @test sp isa Pair
end

# ── ExpandScalarProduct ──────────────────────────────────────────────

@testset "Expand FourVector" begin
    # (p+q)^μ → p^μ + q^μ
    # Ref: ExpandScalarProduct[FV[p+q, mu]]
    pq = Momentum(:p) + Momentum(:q)
    fv = Pair(LorentzIndex(:μ), pq)
    result = expand_scalar_product(fv)
    @test length(result) == 2
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[FV(:p, :μ)] == 1 // 1
    @test coeffs[FV(:q, :μ)] == 1 // 1
end

@testset "Expand ScalarProduct" begin
    # (p+q)·k → p·k + q·k
    # Ref: ExpandScalarProduct[SP[p+q, k]]
    pq = Momentum(:p) + Momentum(:q)
    sp = Pair(Momentum(:k), pq)
    result = expand_scalar_product(sp)
    @test length(result) == 2
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[SP(:k, :p)] == 1 // 1
    @test coeffs[SP(:k, :q)] == 1 // 1
end

@testset "Expand bilinear (both slots)" begin
    # (p+q)·(k+l) → p·k + p·l + q·k + q·l
    pq = Momentum(:p) + Momentum(:q)
    kl = Momentum(:k) + Momentum(:l)
    sp = Pair(pq, kl)
    result = expand_scalar_product(sp)
    @test length(result) == 4
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[SP(:k, :p)] == 1 // 1
    @test coeffs[SP(:l, :p)] == 1 // 1
    @test coeffs[SP(:k, :q)] == 1 // 1
    @test coeffs[SP(:l, :q)] == 1 // 1
end

@testset "Expand with coefficients" begin
    # (2p - q)·k → 2(p·k) - (q·k)
    p2mq = 2 * Momentum(:p) - Momentum(:q)
    sp = Pair(Momentum(:k), p2mq)
    result = expand_scalar_product(sp)
    @test length(result) == 2
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[SP(:k, :p)] == 2 // 1
    @test coeffs[SP(:k, :q)] == -1 // 1
end

@testset "No expansion needed" begin
    # Plain SP (no MomentumSum) → unchanged
    sp = SP(:p, :q)
    result = expand_scalar_product(sp)
    @test length(result) == 1
    @test result[1] == (1 // 1, sp)
end

@testset "Expand with SPContext" begin
    # (p+q)·p with SP[p,p]=m², SP[p,q]=s → m² + s
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :p, :m²)
    ctx = set_sp(ctx, :p, :q, :s)
    pq = Momentum(:p) + Momentum(:q)
    sp = Pair(Momentum(:p), pq)
    result = expand_scalar_product(sp; ctx=ctx)
    @test length(result) == 2
    vals = Dict(r[2] => r[1] for r in result)
    @test vals[:m²] == 1 // 1
    @test vals[:s] == 1 // 1
end

@testset "Expand D-dimensional" begin
    # (p+q)^μ_D → p^μ_D + q^μ_D
    pq = Momentum(:p, DimD()) + Momentum(:q, DimD())
    fv = Pair(LorentzIndex(:μ, DimD()), pq)
    result = expand_scalar_product(fv)
    @test length(result) == 2
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[FVD(:p, :μ)] == 1 // 1
    @test coeffs[FVD(:q, :μ)] == 1 // 1
end
