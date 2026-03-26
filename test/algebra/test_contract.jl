# Tests for Lorentz index contraction
# Translated from FeynCalc Tests/Lorentz/Contract.test and PairContract.test
import Feynfeld: Pair

# ── Metric trace ─────────────────────────────────────────────────────

@testset "Metric trace" begin
    # g^{μ}_{μ} = 4 in 4D
    # Ref: PairContract.test, PairContract[LorentzIndex[i], LorentzIndex[i]] = 4
    @test contract(MT(:μ, :μ)) == 4

    # g^{μ}_{μ} = D in D dimensions
    # Ref: PairContract.test, PairContract[LorentzIndex[i,D], LorentzIndex[i,D]] = D
    @test contract(MTD(:μ, :μ)) == :D

    # g^{μ}_{μ} = D-4 in evanescent dimensions
    @test contract(MTE(:μ, :μ)) == :(D - 4)
end

# ── Metric-metric contraction ────────────────────────────────────────

@testset "Metric-metric contraction" begin
    # g^{μν} g_{νρ} = g^{μ}_{ρ} = δ^μ_ρ
    # Ref: Contract.test, Contract[MT[mu,nu] MT[nu,rho]]
    result = contract(MT(:μ, :ν), MT(:ν, :ρ))
    @test result isa Pair
    @test result == MT(:μ, :ρ)

    # D-dimensional: g^{μν}_D g^D_{νρ} = g^{μρ}_D
    result = contract(MTD(:μ, :ν), MTD(:ν, :ρ))
    @test result == MTD(:μ, :ρ)

    # Chain: g^{μν} g_{νρ} g^{ρσ} = g^{μσ}
    result = contract(MT(:μ, :ν), MT(:ν, :ρ), MT(:ρ, :σ))
    @test result isa Pair
    @test result == MT(:μ, :σ)

    # Full contraction: g^{μν} g_{μν} = 4
    @test contract(MT(:μ, :ν), MT(:μ, :ν)) == 4

    # Full contraction D-dim: g^{μν}_D g^D_{μν} = D
    @test contract(MTD(:μ, :ν), MTD(:μ, :ν)) == :D
end

# ── Metric-vector contraction ────────────────────────────────────────

@testset "Metric-vector contraction" begin
    # g^{μν} p_ν = p^μ
    # Ref: Contract.test
    result = contract(MT(:μ, :ν), FV(:p, :ν))
    @test result isa Pair
    @test result == FV(:p, :μ)

    # D-dim version
    result = contract(MTD(:μ, :ν), FVD(:p, :ν))
    @test result == FVD(:p, :μ)
end

# ── Vector-vector contraction ────────────────────────────────────────

@testset "Vector-vector contraction" begin
    # p^μ q_μ = p·q
    # Ref: Contract[FV[p,mu] FV[q,mu]] = SP[p,q]
    result = contract(FV(:p, :μ), FV(:q, :μ))
    @test result isa Pair
    @test result == SP(:p, :q)

    # D-dim: p^μ_D q^D_μ = SPD[p,q]
    result = contract(FVD(:p, :μ), FVD(:q, :μ))
    @test result == SPD(:p, :q)

    # Self-contraction: p^μ p_μ = p²
    result = contract(FV(:p, :μ), FV(:p, :μ))
    @test result == SP(:p)
end

# ── BMHV mixed-dimension contraction ─────────────────────────────────

@testset "BMHV mixed-dim traces" begin
    # All 9 combinations of MT/MTD/MTE trace contractions
    # Ref: Contract.test BMHV-ID29 through BMHV-ID36
    @test contract(MT(:μ, :ν), MTD(:μ, :ν)) == 4      # 4∩D = 4, trace = 4
    @test contract(MT(:μ, :ν), MTE(:μ, :ν)) == 0       # 4∩(D-4) = 0
    @test contract(MTD(:μ, :ν), MT(:μ, :ν)) == 4       # D∩4 = 4, trace = 4
    @test contract(MTD(:μ, :ν), MTE(:μ, :ν)) == :(D - 4)  # D∩(D-4) = D-4
    @test contract(MTE(:μ, :ν), MT(:μ, :ν)) == 0       # (D-4)∩4 = 0
    @test contract(MTE(:μ, :ν), MTD(:μ, :ν)) == :(D - 4)  # (D-4)∩D = D-4
end

@testset "BMHV metric-metric partial contraction" begin
    # Ref: Contract.test BMHV-ID19 through BMHV-ID27
    @test contract(MT(:μ, :ν), MTD(:ν, :ρ)) == MT(:μ, :ρ)   # 4∩D → 4
    @test contract(MT(:μ, :ν), MTE(:ν, :ρ)) == 0              # 4∩(D-4) → 0
    @test contract(MTD(:μ, :ν), MT(:ν, :ρ)) == MT(:μ, :ρ)    # D∩4 → 4
    @test contract(MTD(:μ, :ν), MTE(:ν, :ρ)) == MTE(:μ, :ρ)  # D∩(D-4) → D-4
    @test contract(MTE(:μ, :ν), MT(:ν, :ρ)) == 0              # (D-4)∩4 → 0
    @test contract(MTE(:μ, :ν), MTD(:ν, :ρ)) == MTE(:μ, :ρ)  # (D-4)∩D → D-4
end

@testset "BMHV metric-vector contraction" begin
    # D-dim metric with 4-dim vector → 4-dim result
    result = contract(MTD(:μ, :ν), FV(:p, :ν))
    @test result == FV(:p, :μ)

    # Evanescent metric with D-dim vector → evanescent
    result = contract(MTE(:μ, :ν), FVD(:p, :ν))
    @test result == FVE(:p, :μ)

    # 4-dim metric with (D-4)-dim vector → vanishes
    @test contract(MT(:μ, :ν), FVE(:p, :ν)) == 0
end

@testset "BMHV vector-vector contraction" begin
    # All non-trivial cross-dim FV contractions
    # Ref: Contract.test BMHV-ID37 through BMHV-ID45
    @test contract(FVD(:p, :μ), FV(:q, :μ)) == SP(:p, :q)     # D∩4 → 4
    @test contract(FV(:p, :μ), FVD(:q, :μ)) == SP(:p, :q)     # 4∩D → 4
    @test contract(FVD(:p, :μ), FVE(:q, :μ)) == SPE(:p, :q)   # D∩(D-4) → D-4
    @test contract(FVE(:p, :μ), FVD(:q, :μ)) == SPE(:p, :q)   # (D-4)∩D → D-4
    @test contract(FV(:p, :μ), FVE(:q, :μ)) == 0               # 4∩(D-4) → 0
    @test contract(FVE(:p, :μ), FV(:q, :μ)) == 0               # (D-4)∩4 → 0
    @test contract(FVE(:p, :μ), FVE(:q, :μ)) == SPE(:p, :q)   # (D-4)∩(D-4)
end

# ── With SPContext ───────────────────────────────────────────────────

@testset "Contract with SPContext" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :p, :m²)
    ctx = set_sp(ctx, :p, :q, :s)

    # p^μ p_μ = m² (from context)
    @test contract(FV(:p, :μ), FV(:p, :μ); ctx=ctx) == :m²

    # p^μ q_μ = s (from context)
    @test contract(FV(:p, :μ), FV(:q, :μ); ctx=ctx) == :s

    # Without context, returns the Pair
    @test contract(FV(:p, :μ), FV(:p, :μ)) == SP(:p)

    # BMHV vanishing takes priority over SP lookup
    @test contract(FV(:p, :μ), FVE(:q, :μ); ctx=ctx) == 0
end

# ── Einstein violation ───────────────────────────────────────────────

@testset "Einstein violation detection" begin
    # Index appearing 3+ times should error
    @test_throws ErrorException contract(
        MT(:μ, :ν), FV(:p, :μ), FV(:q, :μ))
end

# ── No contraction (free indices only) ───────────────────────────────

@testset "No dummy indices" begin
    # Single Pair with distinct indices — nothing to contract
    result = contract(MT(:μ, :ν))
    @test result isa Pair
    @test result == MT(:μ, :ν)
end
