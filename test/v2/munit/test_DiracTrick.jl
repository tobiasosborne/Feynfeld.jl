# MUnit translations: FeynCalc DiracTrick.test
# Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test
# Ground truth: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.9)
# "gamma^mu Gamma^(l) gamma_mu = (-1)^l {(D-2l)Gamma^(l)
#   - 4 sum_{i<j} (-1)^{j-i} Gamma_{ij}^(l) g_{mu_i mu_j}}"
# Convention: all identities verified symbolically (exact AlgSum equality).
#
# Scope: Pure same-dimension tests (GAD-only, GA-only) and GAD-bracket
# with uniform-dimension inner gammas. Tests involving GAE/GSD/GSE or
# BMHV mixed-dimension decomposition are deferred pending those constructors.

using Test
using Feynfeld

@testset "MUnit DiracTrick" begin

    # ---- Helpers ----
    function trick(gammas::Vector{DiracGamma})
        de = DiracExpr([(alg(1), DiracChain(gammas))])
        dirac_trick(de)
    end

    function result_terms(de::DiracExpr)
        [(c, DiracGamma[g for g in ch.elements if g isa DiracGamma]) for (c, ch) in de.terms]
    end

    # Find the (coefficient, chain) pair whose chain matches `target`
    function find_term(terms, target::Vector{DiracGamma})
        for (c, gs) in terms
            gs == target && return c
        end
        nothing
    end

    # ========================================================================
    # OneFreeIndex (n=1): gamma^mu gamma^a gamma_mu = (2-D) gamma^a
    # Ref: MertigBohmDenner1991, Eq. (2.9), l=1
    # ========================================================================
    @testset "OneFreeIndex" begin
        # fcstDiracTrickOneFreeIndex-ID1
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, OneFreeIndex-ID1
        # Input:  GAD[mu, nu, mu]
        # Output: (2 - D)*GAD[nu]
        @testset "ID1: gamma^mu_D gamma^nu_D gamma_{mu,D} = (2-D)gamma^nu_D" begin
            result = trick(DiracGamma[GAD(:mu), GAD(:nu), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 1
            c, gs = terms[1]
            @test gs == DiracGamma[GAD(:nu)]
            @test c == alg(2 - DIM)
        end

        # fcstDiracTrickOneFreeIndex-ID2
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, OneFreeIndex-ID2
        # Input:  GAD[mu].GA[nu].GAD[mu]
        # Output: (2 - D)*GA[nu]
        @testset "ID2: gamma^mu_D gamma^nu_4 gamma_{mu,D} = (2-D)gamma^nu_4" begin
            result = trick(DiracGamma[GAD(:mu), GA(:nu), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 1
            c, gs = terms[1]
            @test gs == DiracGamma[GA(:nu)]
            @test c == alg(2 - DIM)
        end

        # fcstDiracTrickOneFreeIndex-ID4
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, OneFreeIndex-ID4
        # Input:  GA[mu, nu, mu]
        # Output: -2*GA[nu]
        @testset "ID4: gamma^mu_4 gamma^nu_4 gamma_{mu,4} = -2 gamma^nu_4" begin
            result = trick(DiracGamma[GA(:mu), GA(:nu), GA(:mu)])
            terms = result_terms(result)
            @test length(terms) == 1
            c, gs = terms[1]
            @test gs == DiracGamma[GA(:nu)]
            @test c == alg(-2)
        end
    end

    # ========================================================================
    # TwoFreeIndices (n=2): gamma^mu gamma^a gamma^b gamma_mu
    #   = 4 g^{ab} + (D-4) gamma^a gamma^b
    # Ref: MertigBohmDenner1991, Eq. (2.9), l=2
    # ========================================================================
    @testset "TwoFreeIndices" begin
        # fcstDiracTrickTwoFreeIndices-ID1
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, TwoFreeIndices-ID1
        # Input:  GAD[mu, nu, rho, mu]
        # Output: (D-4)*GAD[nu].GAD[rho] + 4*g^{nu,rho}_D
        @testset "ID1: D-dim n=2" begin
            result = trick(DiracGamma[GAD(:mu), GAD(:nu), GAD(:rho), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 2

            chain_fwd = DiracGamma[GAD(:nu), GAD(:rho)]
            chain_empty = DiracGamma[]

            c_fwd = find_term(terms, chain_fwd)
            c_met = find_term(terms, chain_empty)

            @test c_fwd !== nothing
            @test c_met !== nothing

            # Gamma chain coefficient: (D-4)
            @test c_fwd == alg(DIM - 4)
            # Metric term: 4 * g^{nu,rho}_D
            @test c_met == 4//1 * alg(MTD(:nu, :rho))
        end

        # fcstDiracTrickTwoFreeIndices-ID2
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, TwoFreeIndices-ID2
        # Input:  GAD[mu].GA[nu, rho].GAD[mu]
        # Output: (D-4)*GA[nu].GA[rho] + 4*g^{nu,rho}_4
        @testset "ID2: D-dim bracket, 4-dim inner" begin
            result = trick(DiracGamma[GAD(:mu), GA(:nu), GA(:rho), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 2

            chain_fwd = DiracGamma[GA(:nu), GA(:rho)]
            chain_empty = DiracGamma[]

            c_fwd = find_term(terms, chain_fwd)
            c_met = find_term(terms, chain_empty)

            @test c_fwd !== nothing
            @test c_met !== nothing

            @test c_fwd == alg(DIM - 4)
            @test c_met == 4//1 * alg(MT(:nu, :rho))
        end

        # fcstDiracTrickTwoFreeIndices-ID4
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, TwoFreeIndices-ID4
        # Input:  GA[mu, nu, rho, mu]
        # Output: 4*g^{nu,rho}_4
        # Note: (4-4)=0 kills the gamma chain term, only metric survives.
        @testset "ID4: 4-dim n=2 (metric only)" begin
            result = trick(DiracGamma[GA(:mu), GA(:nu), GA(:rho), GA(:mu)])
            terms = result_terms(result)
            @test length(terms) == 1

            c, gs = terms[1]
            @test gs == DiracGamma[]
            @test c == 4//1 * alg(MT(:nu, :rho))
        end
    end

    # ========================================================================
    # ThreeFreeIndices (n=3): gamma^mu gamma^a gamma^b gamma^c gamma_mu
    #   = -(D-4) gamma^a gamma^b gamma^c - 2 gamma^c gamma^b gamma^a
    # Ref: MertigBohmDenner1991, Eq. (2.9), l=3
    # ========================================================================
    @testset "ThreeFreeIndices" begin
        # fcstDiracTrickThreeFreeIndices-ID1
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, ThreeFreeIndices-ID1
        # Input:  GAD[mu, nu, rho, si, mu]
        # Output: -(D-4)*GAD[nu].GAD[rho].GAD[si] - 2*GAD[si].GAD[rho].GAD[nu]
        @testset "ID1: D-dim n=3" begin
            result = trick(DiracGamma[GAD(:mu), GAD(:nu), GAD(:rho), GAD(:si), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 2

            fwd = DiracGamma[GAD(:nu), GAD(:rho), GAD(:si)]
            rev = DiracGamma[GAD(:si), GAD(:rho), GAD(:nu)]

            c_fwd = find_term(terms, fwd)
            c_rev = find_term(terms, rev)

            @test c_fwd !== nothing
            @test c_rev !== nothing

            # Forward: -(D-4) = 4-D
            @test c_fwd == alg(4 - DIM)
            # Reversed: -2
            @test c_rev == alg(-2)
        end

        # fcstDiracTrickThreeFreeIndices-ID2
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, ThreeFreeIndices-ID2
        # Input:  GAD[mu].GA[nu, rho, si].GAD[mu]
        # Output: -(D-4)*GA[nu].GA[rho].GA[si] - 2*GA[si].GA[rho].GA[nu]
        @testset "ID2: D-dim bracket, 4-dim inner" begin
            result = trick(DiracGamma[GAD(:mu), GA(:nu), GA(:rho), GA(:si), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 2

            fwd = DiracGamma[GA(:nu), GA(:rho), GA(:si)]
            rev = DiracGamma[GA(:si), GA(:rho), GA(:nu)]

            c_fwd = find_term(terms, fwd)
            c_rev = find_term(terms, rev)

            @test c_fwd !== nothing
            @test c_rev !== nothing

            @test c_fwd == alg(4 - DIM)
            @test c_rev == alg(-2)
        end

        # fcstDiracTrickThreeFreeIndices-ID4
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, ThreeFreeIndices-ID4
        # Input:  GA[mu, nu, rho, si, mu]
        # Output: -2*GA[si].GA[rho].GA[nu]
        # Note: -(4-4)=0 kills forward term.
        @testset "ID4: 4-dim n=3 (reversed only)" begin
            result = trick(DiracGamma[GA(:mu), GA(:nu), GA(:rho), GA(:si), GA(:mu)])
            terms = result_terms(result)
            @test length(terms) == 1

            c, gs = terms[1]
            @test gs == DiracGamma[GA(:si), GA(:rho), GA(:nu)]
            @test c == alg(-2)
        end
    end

    # ========================================================================
    # FourFreeIndices (n=4): gamma^mu gamma^a gamma^b gamma^c gamma^d gamma_mu
    #   = (D-4) gamma^a gamma^b gamma^c gamma^d
    #   + 2 gamma^c gamma^b gamma^a gamma^d
    #   + 2 gamma^d gamma^a gamma^b gamma^c
    # Ref: MertigBohmDenner1991, Eq. (2.9), l=4
    # ========================================================================
    @testset "FourFreeIndices" begin
        # fcstDiracTrickFourFreeIndices-ID1
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, FourFreeIndices-ID1
        # Input:  GAD[mu, nu, rho, si, tau, mu]
        @testset "ID1: D-dim n=4" begin
            result = trick(DiracGamma[
                GAD(:mu), GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 3

            orig  = DiracGamma[GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau)]
            perm1 = DiracGamma[GAD(:si), GAD(:rho), GAD(:nu), GAD(:tau)]
            perm2 = DiracGamma[GAD(:tau), GAD(:nu), GAD(:rho), GAD(:si)]

            c_orig  = find_term(terms, orig)
            c_perm1 = find_term(terms, perm1)
            c_perm2 = find_term(terms, perm2)

            @test c_orig  !== nothing
            @test c_perm1 !== nothing
            @test c_perm2 !== nothing

            @test c_orig  == alg(DIM - 4)
            @test c_perm1 == alg(2)
            @test c_perm2 == alg(2)
        end

        # fcstDiracTrickFourFreeIndices-ID2
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, FourFreeIndices-ID2
        # Input:  GAD[mu].GA[nu, rho, si, tau].GAD[mu]
        @testset "ID2: D-dim bracket, 4-dim inner" begin
            result = trick(DiracGamma[
                GAD(:mu), GA(:nu), GA(:rho), GA(:si), GA(:tau), GAD(:mu)])
            terms = result_terms(result)
            @test length(terms) == 3

            orig  = DiracGamma[GA(:nu), GA(:rho), GA(:si), GA(:tau)]
            perm1 = DiracGamma[GA(:si), GA(:rho), GA(:nu), GA(:tau)]
            perm2 = DiracGamma[GA(:tau), GA(:nu), GA(:rho), GA(:si)]

            c_orig  = find_term(terms, orig)
            c_perm1 = find_term(terms, perm1)
            c_perm2 = find_term(terms, perm2)

            @test c_orig  !== nothing
            @test c_perm1 !== nothing
            @test c_perm2 !== nothing

            @test c_orig  == alg(DIM - 4)
            @test c_perm1 == alg(2)
            @test c_perm2 == alg(2)
        end

        # fcstDiracTrickFourFreeIndices-ID4
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, FourFreeIndices-ID4
        # Input:  GA[mu, nu, rho, si, tau, mu]
        # Output: 2*perm1 + 2*perm2 (D-4 = 0 kills original)
        @testset "ID4: 4-dim n=4 (two permutations)" begin
            result = trick(DiracGamma[
                GA(:mu), GA(:nu), GA(:rho), GA(:si), GA(:tau), GA(:mu)])
            terms = result_terms(result)
            @test length(terms) == 2

            perm1 = DiracGamma[GA(:si), GA(:rho), GA(:nu), GA(:tau)]
            perm2 = DiracGamma[GA(:tau), GA(:nu), GA(:rho), GA(:si)]

            c_perm1 = find_term(terms, perm1)
            c_perm2 = find_term(terms, perm2)

            @test c_perm1 !== nothing
            @test c_perm2 !== nothing

            @test c_perm1 == alg(2)
            @test c_perm2 == alg(2)
        end
    end

    # ========================================================================
    # FiveFreeIndices (n=5): general formula Eq. (2.9)
    # gamma^mu gamma^{a1}...gamma^{a5} gamma_mu
    #   = -(D-10) gamma^{a1}...gamma^{a5}
    #   + 4 sum_{i<j} (-1)^{j-i+1} g^{ai,aj} [remaining 3 gammas]
    # Ref: MertigBohmDenner1991, Eq. (2.9), l=5
    # ========================================================================
    @testset "FiveFreeIndices" begin
        # fcstDiracTrickFiveFreeIndices-ID1
        # Source: refs/FeynCalc/Tests/Dirac/DiracTrick.test, FiveFreeIndices-ID1
        # Input:  GAD[mu, nu, rho, si, tau, ka, mu]
        # Output: -(D-10)*forward + 4*(10 metric pair terms)
        @testset "ID1: D-dim n=5 general formula" begin
            result = trick(DiracGamma[
                GAD(:mu), GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau), GAD(:ka), GAD(:mu)])
            terms = result_terms(result)

            # 1 forward (5-gamma) + 10 pair (3-gamma) = 11 terms
            @test length(terms) == 11

            # Forward term: -(D-10) = 10-D
            fwd = DiracGamma[GAD(:nu), GAD(:rho), GAD(:si), GAD(:tau), GAD(:ka)]
            c_fwd = find_term(terms, fwd)
            @test c_fwd !== nothing
            @test c_fwd == alg(10 - DIM)

            # Spot-check 3 pair terms from the 10-term sum:

            # Pair (nu,rho): coeff -4*g(nu,rho)_D, chain [si,tau,ka]
            chain_12 = DiracGamma[GAD(:si), GAD(:tau), GAD(:ka)]
            c_12 = find_term(terms, chain_12)
            @test c_12 !== nothing
            @test c_12 == -4//1 * alg(MTD(:nu, :rho))

            # Pair (nu,ka): coeff +4*g(nu,ka)_D, chain [rho,si,tau]
            chain_15 = DiracGamma[GAD(:rho), GAD(:si), GAD(:tau)]
            c_15 = find_term(terms, chain_15)
            @test c_15 !== nothing
            @test c_15 == 4//1 * alg(MTD(:ka, :nu))

            # Pair (si,tau): coeff -4*g(si,tau)_D, chain [nu,rho,ka]
            chain_34 = DiracGamma[GAD(:nu), GAD(:rho), GAD(:ka)]
            c_34 = find_term(terms, chain_34)
            @test c_34 !== nothing
            @test c_34 == -4//1 * alg(MTD(:si, :tau))

            # Count: 1 five-gamma chain + 10 three-gamma chains
            n5 = count(((c, gs),) -> length(gs) == 5, terms)
            n3 = count(((c, gs),) -> length(gs) == 3, terms)
            @test n5 == 1
            @test n3 == 10
        end
    end

end
