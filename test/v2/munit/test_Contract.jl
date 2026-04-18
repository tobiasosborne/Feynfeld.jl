# MUnit translations: FeynCalc Contract.test (Lorentz contraction)
# Source: refs/FeynCalc/Tests/Lorentz/Contract.test
# Ground truth:
#   Metric contraction: Peskin & Schroeder Eq. (A.3), g^{mu}_{mu} = D
#   Levi-Civita: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.21)
#   "eps^{mu nu rho sigma} eps_{mu' nu' rho' sigma'} = -det|g^mu_{mu'} ...|"
# Convention: all identities verified symbolically (exact AlgSum equality).
#
# Scope: metric×metric, metric×FV, FV×FV, Levi-Civita contractions.
# Already in batch1: ID4-7 (metric/FV basics). New here: ID1-3 (eps×eps), ID8, ID10.

using Test
using Feynfeld

@testset "MUnit Contract" begin

    # ---- Helpers ----
    li(s) = LorentzIndex(s, Dim4())   # 4-dim for Levi-Civita
    liD(s) = LorentzIndex(s, DimD())
    mom(s) = Momentum(s)

    # ==== Levi-Civita contractions ====
    # Ref: MertigBohmDenner1991, Eq. (2.21)

    # fcstContractContractionsIn4dims-ID1
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, ID1
    # eps^{ijkl} eps_{ijkl} = -24
    # Full contraction: all 4 indices shared.
    @testset "ID1: eps*eps full contraction = -24" begin
        e1 = alg(Eps(li(:i), li(:j), li(:k), li(:l)))
        e2 = alg(Eps(li(:i), li(:j), li(:k), li(:l)))
        product = e1 * e2
        # Step 1: eps_contract (pair of Eps → determinant of metrics)
        det_expanded = eps_contract(product)
        # Step 2: contract repeated indices
        result = contract(det_expanded)
        @test result == alg(-24)
    end

    # fcstContractContractionsIn4dims-ID2
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, ID2
    # eps^{ijkl} eps_{ijkm} = -6 g^{lm}
    # 3 shared indices (i,j,k), 1 free pair (l vs m)
    @testset "ID2: eps*eps 3 shared = -6 g^{lm}" begin
        e1 = alg(Eps(li(:i), li(:j), li(:k), li(:l)))
        e2 = alg(Eps(li(:i), li(:j), li(:k), li(:m)))
        product = e1 * e2
        det_expanded = eps_contract(product)
        result = contract(det_expanded)
        @test result == -6//1 * alg(pair(li(:l), li(:m)))
    end

    # fcstContractContractionsIn4dims-ID10
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, ID10
    # eps^{a,nu,rho,sigma} eps_{b,nu,rho,sigma} = -6 g^{ab}
    # Same as ID2 with different index labels
    @testset "ID10: eps*eps 3 shared (relabeled)" begin
        e1 = alg(Eps(li(:a), li(:nu), li(:rho), li(:sigma)))
        e2 = alg(Eps(li(:b), li(:nu), li(:rho), li(:sigma)))
        product = e1 * e2
        det_expanded = eps_contract(product)
        result = contract(det_expanded)
        @test result == -6//1 * alg(pair(li(:a), li(:b)))
    end

    # fcstContractContractionsIn4dims-ID3
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, ID3
    # eps^{ijkl} eps_{ijmn} = 2*(g^{kn}g^{lm} - g^{km}g^{ln})
    # 2 shared indices (i,j), 2 free pairs (k,l vs m,n)
    @testset "ID3: eps*eps 2 shared = antisymmetric metric product" begin
        e1 = alg(Eps(li(:i), li(:j), li(:k), li(:l)))
        e2 = alg(Eps(li(:i), li(:j), li(:m), li(:n)))
        product = e1 * e2
        det_expanded = eps_contract(product)
        result = contract(det_expanded)
        expected = 2//1 * (alg(pair(li(:k), li(:n))) * alg(pair(li(:l), li(:m))) -
                           alg(pair(li(:k), li(:m))) * alg(pair(li(:l), li(:n))))
        @test result == expected
    end

    # ==== D-dimensional Levi-Civita contractions ====
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, fcstContractDDims

    # fcstContractDDims-ID2
    # LCD[i1,i2,i3,i4] LCD[i1,i2,i3,i4] = D(1-D)(2-D)(3-D)
    @testset "DDims-ID2: LCD*LCD full = D(1-D)(2-D)(3-D)" begin
        e1 = alg(Eps(liD(:i1), liD(:i2), liD(:i3), liD(:i4)))
        product = e1 * e1
        result = contract(eps_contract(product))
        @test result == alg(DIM * (1 - DIM) * (2 - DIM) * (3 - DIM))
    end

    # fcstContractDDims-ID4
    # LCD[i1,i2,i3,i4] LCD[i1,i2,i3,j4] = (1-D)(2-D)(3-D) g^{i4,j4}_D
    @testset "DDims-ID4: LCD*LCD 3 shared = (1-D)(2-D)(3-D) g_D" begin
        e1 = alg(Eps(liD(:i1), liD(:i2), liD(:i3), liD(:i4)))
        e2 = alg(Eps(liD(:i1), liD(:i2), liD(:i3), liD(:j4)))
        result = contract(eps_contract(e1 * e2))
        coeff = (1 - DIM) * (2 - DIM) * (3 - DIM)
        @test result == alg(coeff) * alg(pair(liD(:i4), liD(:j4)))
    end

    # ==== Additional metric/FV contractions (not in batch1) ====

    # fcstContractContractionsIn4dims-ID8
    # Source: refs/FeynCalc/Tests/Lorentz/Contract.test, ID8
    # (2p)^mu (2p)_mu = 4*p.p
    @testset "ID8: (2p)^mu (2p)_mu = 4*p.p" begin
        # 2p as MomentumSum: 2*Momentum(:p)
        fv_2p = alg(pair(liD(:mu), 2 * mom(:p)))
        product = fv_2p * fv_2p
        # expand_scalar_product first (to resolve MomentumSum)
        expanded = expand_scalar_product(product)
        result = contract(expanded)
        @test result == 4//1 * alg(pair(mom(:p), mom(:p)))
    end

    # ==== Metric contraction sanity (D-dim, complements batch1) ====

    # g^{mu,nu} g_{nu,rho} = g^{mu}_{rho} (D-dim)
    @testset "metric chain: g^{mu,nu} g_{nu,rho} = g^{mu,rho}" begin
        s = alg(pair(liD(:mu), liD(:nu))) * alg(pair(liD(:nu), liD(:rho)))
        result = contract(s)
        @test result == alg(pair(liD(:mu), liD(:rho)))
    end

    # g^{mu,nu} g_{mu,nu} = D (D-dim trace, variant of ID5)
    @testset "D-dim metric trace g^{mu,nu} g_{mu,nu} = D" begin
        s = alg(pair(liD(:mu), liD(:nu))) * alg(pair(liD(:mu), liD(:nu)))
        result = contract(s)
        @test result == alg(DIM)
    end

end
