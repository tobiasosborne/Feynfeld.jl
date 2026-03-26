# Tests for Levi-Civita tensor and EpsContract
# Ref: FeynCalc Tests/Lorentz/EpsContract.test, Eps.test
import Feynfeld: Pair

@testset "Eps construction" begin
    e = LC(:a, :b, :c, :d)
    @test e isa Eps
    @test e.args[1] == LorentzIndex(:a)

    eD = LCD(:a, :b, :c, :d)
    @test eD isa Eps
    @test eD.args[1].dim === DimD()
end

@testset "Eps vanishing" begin
    # Repeated indices → 0 (antisymmetry)
    @test levi_civita(LorentzIndex(:a), LorentzIndex(:a),
                      LorentzIndex(:c), LorentzIndex(:d)) == 0

    # BMHV: (D-4)-dim arg → 0
    @test levi_civita(LorentzIndex(:a, DimDm4()), LorentzIndex(:b),
                      LorentzIndex(:c), LorentzIndex(:d)) == 0

    # Mixed index types with repeat
    @test levi_civita(Momentum(:p), Momentum(:p),
                      LorentzIndex(:c), LorentzIndex(:d)) == 0
end

@testset "Eps with Momentum args" begin
    e = levi_civita(LorentzIndex(:μ), LorentzIndex(:ν), Momentum(:p), Momentum(:q))
    @test e isa Eps
    @test e.args[3] == Momentum(:p)
end

# ── EpsContract: full contraction (4 shared) ─────────────────────────

@testset "EpsContract: 4 shared (4D)" begin
    # ε^{abcd} ε_{abcd} = -24
    # Ref: EpsContract.test ID2, LC[i1,i2,i3,i4]^2
    e = LC(:a, :b, :c, :d)
    @test eps_contract(e, e) == -24
end

@testset "EpsContract: 4 shared (D-dim)" begin
    # ε^{abcd}_D ε^D_{abcd} = -D(D-1)(D-2)(D-3)
    # Ref: EpsContract.test ID5
    e = LCD(:a, :b, :c, :d)
    result = eps_contract(e, e)
    @test result isa Expr
end

# ── EpsContract: 3 shared ────────────────────────────────────────────

@testset "EpsContract: 3 shared (4D)" begin
    # ε^{abcμ} ε_{abcν} = -6 g^{μν}
    # Ref: EpsContract.test ID17
    e1 = LC(:a, :b, :c, :μ)
    e2 = LC(:a, :b, :c, :ν)
    result = eps_contract(e1, e2)
    @test result == (-6, MT(:μ, :ν))

    # Order shouldn't matter for the shared indices
    e1b = LC(:c, :a, :b, :μ)
    e2b = LC(:a, :c, :b, :ν)
    result2 = eps_contract(e1b, e2b)
    @test result2 == (-6, MT(:μ, :ν))
end

@testset "EpsContract: 3 shared BMHV vanishing" begin
    # 4D eps with (D-4) free index → vanishes via pair()
    e1 = LC(:a, :b, :c, :μ)
    e2_args = (LorentzIndex(:a), LorentzIndex(:b), LorentzIndex(:c),
               LorentzIndex(:ν, DimDm4()))
    e2 = Eps(e2_args)
    @test eps_contract(e1, e2) == 0
end

# ── EpsContract: 2 shared ────────────────────────────────────────────

@testset "EpsContract: 2 shared (4D)" begin
    # ε^{abμν} ε_{abρσ} = -2(g^{μρ}g^{νσ} - g^{μσ}g^{νρ})
    # Ref: EpsContract.test ID16
    e1 = LC(:a, :b, :μ, :ν)
    e2 = LC(:a, :b, :ρ, :σ)
    result = eps_contract(e1, e2)
    @test result isa Vector
    @test length(result) == 2
    # Check coefficients are -2 and +2
    coeffs = sort([r[1] for r in result])
    @test coeffs == [-2, 2]
end

