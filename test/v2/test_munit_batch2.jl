# MUnit test translation batch 2: DiracTrick n≥3
# Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test
# Ground truth: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9)
# Translated per protocol: CLAUDE.md §MUnit translation protocol.

using Test
@isdefined(FeynfeldX) || include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "MUnit batch 2: DiracTrick n≥3" begin

    # Helper: build DiracExpr from a single chain and apply dirac_trick
    function trick(gammas::Vector{DiracGamma})
        de = DiracExpr([(alg(1), DiracChain(gammas))])
        dirac_trick(de)
    end

    # Helper: extract (coeff, chain_gammas) pairs from DiracExpr result
    function result_terms(de::DiracExpr)
        [(c, DiracGamma[g for g in ch.elements if g isa DiracGamma]) for (c, ch) in de.terms]
    end

    # ==== ThreeFreeIndices ====

    @testset "DiracTrick n=3 (ThreeFreeIndices)" begin
        # fcstDiracTrickThreeFreeIndices-ID1
        # Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9), l=3
        # "γ^μ γ^ν γ^ρ γ^σ γ_μ = -(D-4) γ^ν γ^ρ γ^σ - 2 γ^σ γ^ρ γ^ν"
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, ThreeFreeIndices-ID1
        # Input:  DiracTrick[GAD[mu, nu, rho, si, mu]]
        # Output: -((-4+D) GAD[nu].GAD[rho].GAD[si]) - 2 GAD[si].GAD[rho].GAD[nu]

        result = trick(DiracGamma[GAD(:mu), GAD(:nu), GAD(:rho), GAD(:si), GAD(:mu)])
        terms = result_terms(result)

        @test length(terms) == 2

        # Identify terms by chain content
        fwd = DiracGamma[GAD(:nu), GAD(:rho), GAD(:si)]   # forward order
        rev = DiracGamma[GAD(:si), GAD(:rho), GAD(:nu)]   # reversed order

        fwd_coeff = nothing
        rev_coeff = nothing
        for (c, gs) in terms
            if gs == fwd
                fwd_coeff = c
            elseif gs == rev
                rev_coeff = c
            end
        end

        @test fwd_coeff !== nothing
        @test rev_coeff !== nothing

        # Forward term: -(D-4) = (4-D)
        if fwd_coeff !== nothing
            fk = FactorKey()
            @test haskey(fwd_coeff.terms, fk)
            @test fwd_coeff.terms[fk] == (4 - DIM)  # -(D-4) = 4-D
        end

        # Reversed term: -2
        if rev_coeff !== nothing
            fk = FactorKey()
            @test haskey(rev_coeff.terms, fk)
            @test rev_coeff.terms[fk] == -2//1
        end

        # Numerical check at D=4: should give -2 γ^σ γ^ρ γ^ν (forward term vanishes)
        @test length(terms) == 2  # structural check
    end

    # ==== FourFreeIndices ====

    @testset "DiracTrick n=4 (FourFreeIndices)" begin
        # fcstDiracTrickFourFreeIndices-ID1
        # Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9), l=4
        # "γ^μ γ^ν γ^ρ γ^σ γ^τ γ_μ = (D-4) γ^ν γ^ρ γ^σ γ^τ
        #   + 2 γ^σ γ^ρ γ^ν γ^τ + 2 γ^τ γ^ν γ^ρ γ^σ"
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, FourFreeIndices-ID1
        # Input:  DiracTrick[GAD[mu, nu, rho, si, tau, mu]]
        # Output: ((-4+D) GAD[nu].GAD[rho].GAD[si].GAD[tau])
        #         + 2 GAD[si].GAD[rho].GAD[nu].GAD[tau]
        #         + 2 GAD[tau].GAD[nu].GAD[rho].GAD[si]

        result = trick(DiracGamma[
            GAD(:mu), GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau), GAD(:mu)])
        terms = result_terms(result)

        @test length(terms) == 3

        orig  = DiracGamma[GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau)]
        perm1 = DiracGamma[GAD(:si), GAD(:rho), GAD(:nu), GAD(:tau)]
        perm2 = DiracGamma[GAD(:tau), GAD(:nu), GAD(:rho), GAD(:si)]

        orig_c = nothing; perm1_c = nothing; perm2_c = nothing
        for (c, gs) in terms
            if gs == orig;  orig_c  = c
            elseif gs == perm1; perm1_c = c
            elseif gs == perm2; perm2_c = c
            end
        end

        @test orig_c  !== nothing
        @test perm1_c !== nothing
        @test perm2_c !== nothing

        # Original order: (D-4)
        if orig_c !== nothing
            fk = FactorKey()
            @test haskey(orig_c.terms, fk)
            @test orig_c.terms[fk] == (DIM - 4)
        end

        # Permutation terms: coefficient 2
        if perm1_c !== nothing
            fk = FactorKey()
            @test haskey(perm1_c.terms, fk)
            @test perm1_c.terms[fk] == 2//1
        end
        if perm2_c !== nothing
            fk = FactorKey()
            @test haskey(perm2_c.terms, fk)
            @test perm2_c.terms[fk] == 2//1
        end
    end
end
