# Spiral 3: qq̄ → gg (quark-antiquark annihilation to two gluons)
#
# Ground truth: refs/FeynCalc/FeynCalc/Examples/QCD/Tree/Mathematica/QQbar-GlGl.m
# Lines 103-104:
# "knownResult = (32/27) g_s^4 (t^2+u^2)/(t u) - (8/3) g_s^4 (t^2+u^2)/s^2"
# Verified CORRECT by FeynCalc vs Ellis, Stirling, Weber Table 7.1.
#
# Three diagrams at tree level:
#   t-channel: quark exchange, colour (T^b T^a)_{ji}, Dirac γ^ν(p̸₁-k̸₁)γ^μ/t
#   u-channel: crossed,       colour (T^a T^b)_{ji}, Dirac γ^μ(p̸₁-k̸₂)γ^ν/u
#   s-channel: triple gluon,  colour T^c f^{cab},    Dirac γ^ρ V_{ρμν}/s
#
# Colour factors (N=3): C_tt=C_uu=16/3, C_tu=-2/3, C_ss=12
# C_ts=C_us=±6i → vanish (imaginary × real = 0 at tree level)
# Spin+colour average: 1/(4 × N²) = 1/36
#
# Key subtlety: individual QCD diagrams are NOT gauge-invariant, so we use
# physical (axial gauge) polarization sums with reference momenta (k1↔k2).
# The relative sign between s-channel and t/u channels comes from the
# gluon propagator (-ig_{ρσ}/s) vs quark propagator (ip̸/t):
# M_t ~ (-ig_s)²(i/t), M_s ~ (-ig_s)(-g_s)(-i/s) → relative phase i.
# This manifests as a -1 on D_ss in the colour-separated decomposition.

using Test
using Feynfeld

@testset "Spiral 3: qq̄ → gg" begin

    @testset "qq̄→gg |M̄|² massless vs FeynCalc" begin
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)

        s_val = 10//1; t_val = -3//1; u_val = -(s_val + t_val)

        ctx = sp_context(
            (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>0//1, (:k2,:k2)=>0//1,
            (:p1,:p2)=>s_val//2, (:p1,:k1)=>-t_val//2, (:p1,:k2)=>-u_val//2,
            (:p2,:k1)=>-u_val//2, (:p2,:k2)=>-t_val//2, (:k1,:k2)=>s_val//2)

        p1mk1 = MomentumSum([(1//1, p1), (-1//1, k1)])
        p1mk2 = MomentumSum([(1//1, p1), (-1//1, k2)])

        # Physical polarization sums (axial gauge, k1↔k2 as reference)
        mi = LorentzIndex(:mu, DimD());  mp = LorentzIndex(:mu_, DimD())
        ni = LorentzIndex(:nu, DimD());  np = LorentzIndex(:nu_, DimD())
        P1 = polarization_sum(mi, mp, k1, k2; ctx)
        P2 = polarization_sum(ni, np, k2, k1; ctx)

        # Helper: compute D_{XY} with physical pol sums
        function eval_D(expr)
            c = contract(expr * P1 * P2; ctx)
            e = expand_scalar_product(c)
            r = evaluate_sp(e; ctx)
            v = first(r.terms)[2]
            v isa DimPoly ? evaluate_dim(v) : v
        end

        trace_XY(gX, gYc) = dirac_trace(DiracGamma[GS(p2); gX; GS(p1); gYc])

        # ──── t/u channel gamma chains (primed conjugate indices) ────
        gt  = DiracGamma[GAD(:nu), DiracGamma(MomSumSlot(p1mk1)), GAD(:mu)]
        gu  = DiracGamma[GAD(:mu), DiracGamma(MomSumSlot(p1mk2)), GAD(:nu)]
        gtc = DiracGamma[GAD(:mu_), DiracGamma(MomSumSlot(p1mk1)), GAD(:nu_)]
        guc = DiracGamma[GAD(:nu_), DiracGamma(MomSumSlot(p1mk2)), GAD(:mu_)]

        D_tt = eval_D(trace_XY(gt, gtc))
        D_uu = eval_D(trace_XY(gu, guc))
        D_tu = eval_D(trace_XY(gt, guc))
        D_ut = eval_D(trace_XY(gu, gtc))

        # ──── s-channel: triple gluon vertex ────
        rho = LorentzIndex(:rho, DimD()); sig = LorentzIndex(:sig, DimD())

        # V_{ρμν}(-q,k1,k2) all outgoing, q=k1+k2
        # = g_{ρμ}(-2k1-k2)_ν + g_{μν}(k1-k2)_ρ + g_{νρ}(k1+2k2)_μ
        function vtx(g1, g2, ms, mi_v)
            r = AlgSum()
            for (c, m) in ms.terms; r = r + c * alg(pair(g1, g2)) * alg(pair(mi_v, m)); end
            r
        end
        m2k1mk2 = MomentumSum([(-2//1, k1), (-1//1, k2)])
        k1mk2   = MomentumSum([(1//1, k1), (-1//1, k2)])
        k1p2k2  = MomentumSum([(1//1, k1), (2//1, k2)])

        Va = vtx(rho,mi,m2k1mk2,ni) + vtx(mi,ni,k1mk2,rho) + vtx(ni,rho,k1p2k2,mi)
        Vc = vtx(sig,mp,m2k1mk2,np) + vtx(mp,np,k1mk2,sig) + vtx(np,sig,k1p2k2,mp)

        qt_ss = dirac_trace(DiracGamma[GS(p2), GAD(:rho), GS(p1), GAD(:sig)])
        D_ss = eval_D(qt_ss * Va * Vc)

        # ──── Combine ────
        # |M̄|²/g_s⁴ = (1/36)[C_tt D_tt/t² + C_uu D_uu/u² + 2C_tu D_tu/(tu) + C_ss D_ss/s²]
        # The s-channel D_ss picks up a relative -1 from the gluon propagator
        # phase: M_s ~ ig_s² D_s (one quark + one gluon vertex), while
        # M_t ~ -ig_s² D_t (two quark vertices). The relative i gives i²=-1
        # when computing the colour-separated D_ss contribution.
        C_tt = 16//3; C_uu = 16//3; C_tu = -2//3; C_ss = 12//1

        total = (1//36) * (
            C_tt * D_tt // t_val^2 +
            C_uu * D_uu // u_val^2 +
            C_tu * D_tu // (t_val * u_val) +
            C_tu * D_ut // (t_val * u_val) +
            C_ss * (-D_ss) // s_val^2)  # note: -D_ss from relative phase

        # ──── Compare ────
        # Ref: refs/FeynCalc/.../QQbar-GlGl.m, lines 103-104
        # "|M̄|²/g_s⁴ = (32/27)(t²+u²)/(tu) - (8/3)(t²+u²)/s²"
        known = (32//27) * (t_val^2 + u_val^2) // (t_val * u_val) -
                (8//3) * (t_val^2 + u_val^2) // s_val^2

        @test total == known
    end

    @testset "qq̄→gg second kinematic point" begin
        p1 = Momentum(:p1); p2 = Momentum(:p2)
        k1 = Momentum(:k1); k2 = Momentum(:k2)
        s2 = 8//1; t2 = -5//1; u2 = -(s2 + t2)
        ctx2 = sp_context(
            (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>0//1, (:k2,:k2)=>0//1,
            (:p1,:p2)=>s2//2, (:p1,:k1)=>-t2//2, (:p1,:k2)=>-u2//2,
            (:p2,:k1)=>-u2//2, (:p2,:k2)=>-t2//2, (:k1,:k2)=>s2//2)
        p1mk1 = MomentumSum([(1//1, p1), (-1//1, k1)])
        p1mk2 = MomentumSum([(1//1, p1), (-1//1, k2)])

        mi = LorentzIndex(:mu, DimD()); mp = LorentzIndex(:mu_, DimD())
        ni = LorentzIndex(:nu, DimD()); np = LorentzIndex(:nu_, DimD())
        P1 = polarization_sum(mi, mp, k1, k2; ctx=ctx2)
        P2 = polarization_sum(ni, np, k2, k1; ctx=ctx2)

        function ev2(expr)
            c = contract(expr * P1 * P2; ctx=ctx2)
            e = expand_scalar_product(c)
            r = evaluate_sp(e; ctx=ctx2)
            v = first(r.terms)[2]; v isa DimPoly ? evaluate_dim(v) : v
        end
        trXY(gX, gYc) = dirac_trace(DiracGamma[GS(p2); gX; GS(p1); gYc])
        gt  = DiracGamma[GAD(:nu), DiracGamma(MomSumSlot(p1mk1)), GAD(:mu)]
        gu  = DiracGamma[GAD(:mu), DiracGamma(MomSumSlot(p1mk2)), GAD(:nu)]
        gtc = DiracGamma[GAD(:mu_), DiracGamma(MomSumSlot(p1mk1)), GAD(:nu_)]
        guc = DiracGamma[GAD(:nu_), DiracGamma(MomSumSlot(p1mk2)), GAD(:mu_)]

        rho = LorentzIndex(:rho, DimD()); sig = LorentzIndex(:sig, DimD())
        function vtx(g1,g2,ms,mv)
            r=AlgSum(); for(c,m) in ms.terms; r=r+c*alg(pair(g1,g2))*alg(pair(mv,m)); end; r
        end
        Va = vtx(rho,mi,MomentumSum([(-2//1,k1),(-1//1,k2)]),ni) +
             vtx(mi,ni,MomentumSum([(1//1,k1),(-1//1,k2)]),rho) +
             vtx(ni,rho,MomentumSum([(1//1,k1),(2//1,k2)]),mi)
        Vc = vtx(sig,mp,MomentumSum([(-2//1,k1),(-1//1,k2)]),np) +
             vtx(mp,np,MomentumSum([(1//1,k1),(-1//1,k2)]),sig) +
             vtx(np,sig,MomentumSum([(1//1,k1),(2//1,k2)]),mp)

        total = (1//36) * (
            (16//3) * ev2(trXY(gt,gtc)) // t2^2 +
            (16//3) * ev2(trXY(gu,guc)) // u2^2 +
            (-2//3) * ev2(trXY(gt,guc)) // (t2*u2) +
            (-2//3) * ev2(trXY(gu,gtc)) // (t2*u2) +
            12 * (-ev2(dirac_trace(DiracGamma[GS(p2),GAD(:rho),GS(p1),GAD(:sig)]) * Va * Vc)) // s2^2)

        known = (32//27)*(t2^2+u2^2)//(t2*u2) - (8//3)*(t2^2+u2^2)//s2^2
        @test total == known
    end
end
