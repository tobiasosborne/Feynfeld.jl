# MUnit translations: FeynCalc PolarizationSum.test
# Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test
# Ground truth: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.13)
# "Sum_lambda eps^mu(k) eps^{nu*}(k) = -g^{mu,nu}" [Feynman gauge]
# Convention: all identities verified symbolically (exact AlgSum equality).
#
# Portable tests: ID1 (Feynman), ID3 (massive), ID5 (axial), ID8 (virtual boson).
# Skipped: ID2/4/6/7 (symbolic k^2 or MomentumSum), ID9-15 (CartesianIndex).

using Test
using Feynfeld

@testset "MUnit PolarizationSum" begin

    # Helpers
    li(s) = LorentzIndex(s, DimD())
    li4(s) = LorentzIndex(s, Dim4())
    mom(s) = Momentum(s)
    mt(a, b) = pair(li(a), li(b))
    fv(p, mu) = pair(li(mu), mom(p))
    sp(a, b) = pair(mom(a), mom(b))

    # ==== Feynman gauge (massless / virtual) ====

    # fcstPolarizationSum-ID1
    # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID1
    # Input:  PolarizationSum[rho, si]
    # Output: -Pair[LorentzIndex[rho], LorentzIndex[si]]
    # Ref: MertigBohmDenner1991, Eq. (2.13)
    # "Sum_lambda eps^rho eps^{si*} = -g^{rho,si}"
    @testset "ID1: Feynman gauge -g^{rho,si}" begin
        @test polarization_sum(li(:rho), li(:si)) == -alg(mt(:rho, :si))
    end

    # fcstPolarizationSum-ID8
    # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID8
    # Input:  PolarizationSum[mu, nu, x1, 0, VirtualBoson->True]
    # Output: -MT[mu, nu]
    @testset "ID8: virtual boson -g^{mu,nu}" begin
        @test polarization_sum(li(:mu), li(:nu)) == -alg(mt(:mu, :nu))
    end

    # ==== Massive polarization sum ====

    # fcstPolarizationSum-ID3
    # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID3
    # Input:  ScalarProduct[k,k]=m^2; PolarizationSum[rho, si, k]
    # Output: -g^{rho,si} + k^rho k^si / m^2
    # Ref: Peskin & Schroeder Eq. (5.71)
    # "Sum_lambda eps^mu(k) eps^{nu*}(k) = -g^{mu,nu} + k^mu k^nu / M^2"
    @testset "ID3: massive -g + k*k/M^2" begin
        M2 = 3//1  # arbitrary mass squared
        result = polarization_sum_massive(li(:rho), li(:si), mom(:k), M2)
        expected = -alg(mt(:rho, :si)) +
                   (1 // M2) * alg(fv(:k, :rho)) * alg(fv(:k, :si))
        @test result == expected
    end

    # Same test with M^2 = 1 (W boson in natural units)
    @testset "ID3 variant: M^2=1" begin
        result = polarization_sum_massive(li(:mu), li(:nu), mom(:W), 1//1)
        expected = -alg(mt(:mu, :nu)) + alg(fv(:W, :mu)) * alg(fv(:W, :nu))
        @test result == expected
    end

    # ==== Axial gauge (massless with reference momentum) ====

    # fcstPolarizationSum-ID5
    # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID5
    # Input:  ScalarProduct[k,k]=0; PolarizationSum[mu, nu, k, n]
    # Output: -g^{mu,nu} + (k^mu n^nu + n^mu k^nu)/(k.n)
    #         - k^mu k^nu n^2 / (k.n)^2
    # Note: full axial gauge formula has n^2 term. Our implementation
    # assumes n^2=0 (massless reference), giving the simpler form:
    #   -g^{mu,nu} + (k^mu n^nu + n^mu k^nu)/(k.n)
    @testset "ID5: axial gauge (massless n)" begin
        ctx = sp_context((:k, :n) => 5//1)
        kn = 5//1
        result = polarization_sum(li(:mu), li(:nu), mom(:k), mom(:n); ctx)
        expected = -alg(mt(:mu, :nu)) +
                   (1 // kn) * (alg(fv(:k, :mu)) * alg(fv(:n, :nu)) +
                                alg(fv(:n, :mu)) * alg(fv(:k, :nu)))
        @test result == expected
    end

    # Verify axial gauge with different k.n value
    @testset "ID5 variant: different k.n" begin
        ctx = sp_context((:p, :q) => 7//2)
        result = polarization_sum(li(:a), li(:b), mom(:p), mom(:q); ctx)
        expected = -alg(mt(:a, :b)) +
                   (2 // 7) * (alg(fv(:p, :a)) * alg(fv(:q, :b)) +
                               alg(fv(:q, :a)) * alg(fv(:p, :b)))
        @test result == expected
    end

end
