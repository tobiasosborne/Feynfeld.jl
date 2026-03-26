# Phase 1a validation: translated FeynCalc MUnit tests
# Source: refs/FeynCalc/Tests/Lorentz/
# Only tests NOT already covered in test_{pair,contract,expand_sp,eps}.jl
import Feynfeld: Pair

# ── PairContract.test: BMHV projection via pair() ────────────────────

@testset "MUnit PairContract: Momentum BMHV" begin
    # fcstPairContract-ID9: Pair[Momentum[p], Momentum[p]] = SP[p,p]
    @test pair(Momentum(:p), Momentum(:p)) == SP(:p, :p)

    # fcstPairContract-ID10: Pair[Momentum[p,D], Momentum[p,D]] = SPD[p,p]
    @test pair(Momentum(:p, DimD()), Momentum(:p, DimD())) == SPD(:p, :p)

    # fcstPairContract-ID12: Momentum D meets 4 → projects to 4
    @test pair(Momentum(:p, DimD()), Momentum(:p, Dim4())) == SP(:p, :p)

    # fcstPairContract-ID13: Momentum (D-4) meets 4 → vanishes
    @test pair(Momentum(:p, DimDm4()), Momentum(:p, Dim4())) == 0

    # fcstPairContract-ID14: Momentum (D-4) meets D → projects to (D-4)
    @test pair(Momentum(:p, DimDm4()), Momentum(:p, DimD())) == SPE(:p, :p)
end

@testset "MUnit PairContract: zero momentum" begin
    # fcstPairContract2-ID1: PairContract[0, x] = 0
    @test pair(MomentumSum(Tuple{Rational{Int},Symbol}[], Dim4()),
              LorentzIndex(:x)) == 0
end

@testset "MUnit PairContract3: SPContext lookup via contract" begin
    # fcstPairContract3-ID1: SP[p,p] = (a+b)^2
    ctx = set_sp(SPContext(), :p, :p, :((a + b)^2))
    @test contract(FV(:p, :mu), FV(:p, :mu); ctx=ctx) == :((a + b)^2)
end

# ── Contract.test: Levi-Civita via eps_contract ──────────────────────

@testset "MUnit Contract 4D: eps_contract patterns" begin
    # fcstContractContractionsIn4dims-ID1
    @test eps_contract(LC(:i, :j, :k, :l), LC(:i, :j, :k, :l)) == -24

    # fcstContractContractionsIn4dims-ID2
    result = eps_contract(LC(:i, :j, :k, :l), LC(:i, :j, :k, :m))
    @test result == (-6, MT(:l, :m))

    # fcstContractContractionsIn4dims-ID3: 2 shared
    result = eps_contract(LC(:i, :j, :k, :l), LC(:i, :j, :m, :n))
    @test length(result) == 2
    term_sets = [(t[1], Set(t[2])) for t in result]
    @test (2, Set([MT(:k, :n), MT(:l, :m)])) in term_sets
    @test (-2, Set([MT(:k, :m), MT(:l, :n)])) in term_sets

    # fcstContractContractionsIn4dims-ID10
    result = eps_contract(LC(:a, :nu, :rho, :sigma), LC(:b, :nu, :rho, :sigma))
    @test result == (-6, MT(:a, :b))
end

# ── Contract.test BMHV: already covered by test_contract.jl ─────────
# (IDs 10-45 fully covered there)

# ── ExpandScalarProduct.test ─────────────────────────────────────────

@testset "MUnit ExpandSP: nested coefficients" begin
    # fcstExpandScalarProduct-ID4: a + 2*(b + 3c) → a + 2b + 6c
    b3c = Momentum(:b) + 3 * Momentum(:c)
    lhs = Momentum(:a) + 2 * b3c
    fv = Pair(LorentzIndex(:mu), lhs)
    result = expand_scalar_product(fv)
    @test length(result) == 3
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[FV(:a, :mu)] == 1 // 1
    @test coeffs[FV(:b, :mu)] == 2 // 1
    @test coeffs[FV(:c, :mu)] == 6 // 1
end

@testset "MUnit ExpandSP: 3x3 bilinear" begin
    # fcstExpandScalarProduct-ID5: (a+2b+6c)·(e+5f+10g) → 9 terms
    lhs = Momentum(:a) + 2 * (Momentum(:b) + 3 * Momentum(:c))
    rhs = Momentum(:e) + 5 * (Momentum(:f) + 2 * Momentum(:g))
    sp = Pair(lhs, rhs)
    result = expand_scalar_product(sp)
    @test length(result) == 9
    coeffs = Dict(r[2] => r[1] for r in result)
    @test coeffs[SP(:a, :e)] == 1 // 1
    @test coeffs[SP(:a, :f)] == 5 // 1
    @test coeffs[SP(:a, :g)] == 10 // 1
    @test coeffs[SP(:b, :e)] == 2 // 1
    @test coeffs[SP(:b, :f)] == 10 // 1
    @test coeffs[SP(:b, :g)] == 20 // 1
    @test coeffs[SP(:c, :e)] == 6 // 1
    @test coeffs[SP(:c, :f)] == 30 // 1
    @test coeffs[SP(:c, :g)] == 60 // 1
end

@testset "MUnit ExpandSP: two-term SP addition" begin
    # fcstExpandScalarProduct-ID2: SP[a+b, c+d] + SP[a+l, c+d]
    ab = Momentum(:a) + Momentum(:b)
    cd = Momentum(:c) + Momentum(:d)
    al = Momentum(:a) + Momentum(:l)
    r1 = expand_scalar_product(Pair(ab, cd))
    r2 = expand_scalar_product(Pair(al, cd))
    combined = Dict{Any,Rational{Int}}()
    for (c, t) in vcat(r1, r2)
        combined[t] = get(combined, t, 0 // 1) + c
    end
    @test combined[SP(:a, :c)] == 2 // 1
    @test combined[SP(:a, :d)] == 2 // 1
    @test combined[SP(:b, :c)] == 1 // 1
    @test combined[SP(:b, :d)] == 1 // 1
    @test combined[SP(:c, :l)] == 1 // 1
    @test combined[SP(:d, :l)] == 1 // 1
end

@testset "MUnit ExpandSP: SPContext via MomentumSum" begin
    # fcstExpandScalarProduct-ID12: SP[p,p] = (a+b)^2
    ctx = set_sp(SPContext(), :p, :p, :((a + b)^2))
    p_sum = MomentumSum([(1 // 1, :p)], Dim4())
    sp = Pair(p_sum, Momentum(:p))
    result = expand_scalar_product(sp; ctx=ctx)
    @test length(result) == 1
    @test result[1] == (1 // 1, :((a + b)^2))
end

# ── ScalarProduct.test ───────────────────────────────────────────────

@testset "MUnit ScalarProduct: construction identity" begin
    # fcstScalarProduct-ID1 through ID6
    @test SP(:a, :a) == Pair(Momentum(:a), Momentum(:a))
    @test SPD(:a, :a) == Pair(Momentum(:a, DimD()), Momentum(:a, DimD()))
    @test SPE(:a, :a) == Pair(Momentum(:a, DimDm4()), Momentum(:a, DimDm4()))
    @test SP(:a, :b) == Pair(Momentum(:a), Momentum(:b))
    @test SPD(:a, :b) == Pair(Momentum(:a, DimD()), Momentum(:b, DimD()))
    @test SPE(:a, :b) == Pair(Momentum(:a, DimDm4()), Momentum(:b, DimDm4()))
end

@testset "MUnit ScalarProduct: symmetry" begin
    @test SP(:a, :b) == SP(:b, :a)
    @test SPD(:a, :b) == SPD(:b, :a)
    @test SPE(:a, :b) == SPE(:b, :a)
end

@testset "MUnit SPContext: contract integration" begin
    ctx = set_sp(SPContext(), :a, :b, :abval)
    @test contract(FV(:a, :mu), FV(:b, :mu); ctx=ctx) == :abval
    ctx2 = set_sp(SPContext(), :a, :a, :(m^2))
    @test contract(FV(:a, :mu), FV(:a, :mu); ctx=ctx2) == :(m^2)
end

@testset "MUnit SPContext: expand integration" begin
    ctx = set_sp(SPContext(), :a, :b, :abval)
    ak = Momentum(:a) + Momentum(:k)
    sp = Pair(ak, Momentum(:b))
    result = expand_scalar_product(sp; ctx=ctx)
    @test length(result) == 2
    vals = Dict(r[2] => r[1] for r in result)
    @test vals[:abval] == 1 // 1
    @test vals[SP(:b, :k)] == 1 // 1
end
