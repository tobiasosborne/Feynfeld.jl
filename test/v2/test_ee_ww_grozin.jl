# ee→WW multi-channel validation against Grozin formula.
#
# Ground truth:
#   Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.9)-(11.10)
#   Ref: refs/FeynCalc/.../AnelEl-WW.m, lines 80, 95-98, 193-198
#   Ref: src/v2/ew_cross_section.jl (Grozin formula, verified vs Denner Tab. 11.4)
#
# Coupling constants NOT in pipeline output (must multiply explicitly):
#   Ref: Denner1993, Eq. (11.9); PDG2024, Table 10.3
#   c_sγ = e² = 4πα              (eeγ × WWγ)
#   c_sZ = e²/(2sin²θ_W)         (eeZ × WWZ, g_V/g_A already in vertex_structure)
#   c_tν = e²/(2sin²θ_W)         (eνW × eνW, P_L already in vertex_structure)

using Test
@isdefined(FeynfeldX) || include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX
using QuadGK: quadgk

@testset "ee→WW Grozin validation" begin

# ── Shared setup ──
model = ew_model(); rules = feynman_rules(model)
p1 = Momentum(:p1); p2 = Momentum(:p2)
k1 = Momentum(:k1); k2 = Momentum(:k2)
incoming = [ExternalLeg(:e, p1, true, false), ExternalLeg(:e, p2, true, true)]
outgoing = [ExternalLeg(:W, k1, false, false), ExternalLeg(:W, k2, false, true)]
channels = tree_channels(model, rules, incoming, outgoing)

amp_sg = build_amplitude(first(c for c in channels if c.channel==:s && c.exchanged==:gamma), rules, model)
amp_sZ = build_amplitude(first(c for c in channels if c.channel==:s && c.exchanged==:Z), rules, model)
amp_t  = build_amplitude(first(c for c in channels if c.channel==:t), rules, model)

# Coupling constants — Ref: Denner1993, Eq. (11.9)
const_e2   = 4π * EW_ALPHA
const_c_sg = const_e2
const_c_sZ = const_e2 / (2 * EW_SIN2_W)
const_c_t  = const_e2 / (2 * EW_SIN2_W)
const_MZ2  = 1.0 / (1.0 - EW_SIN2_W)  # M_Z² in M_W units
const_MW   = EW_M_W
const_PB   = 3.894e8  # GeV⁻² → pb

# Shared integer kinematic point (M_W units)
s_val = 8//1; t_val = -1//1; u_val = 4//1 - s_val - t_val
ctx_int = sp_context(
    (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>1//1, (:k2,:k2)=>1//1,
    (:p1,:p2)=>s_val//2, (:p1,:k1)=>-t_val//2, (:p1,:k2)=>-u_val//2,
    (:p2,:k1)=>-u_val//2, (:p2,:k2)=>-t_val//2, (:k1,:k2)=>(s_val-2)//2)

# Helper: extract scalar from evaluated AlgSum
_scalar(r) = let c = first(r.terms)[2]; c isa DimPoly ? evaluate_dim(c) : c end

# ── Helper: fully contract gauge exchange expression ──
function gauge_contract(chain_de::DiracExpr, vtx::AlgSum)
    tr = FeynfeldX._single_line_trace(chain_de)
    neg_q = MomentumSum([(-1//1, p1), (-1//1, p2)])
    vtx_c = triple_gauge_vertex(LorentzIndex(:rho_s_, DimD()),
        LorentzIndex(:mu_k1_, DimD()), LorentzIndex(:mu_k2_, DimD()), neg_q, k1, k2)
    mi = LorentzIndex(:mu_k1, DimD()); mp = LorentzIndex(:mu_k1_, DimD())
    ni = LorentzIndex(:mu_k2, DimD()); np = LorentzIndex(:mu_k2_, DimD())
    P1 = polarization_sum_massive(mi, mp, k1, 1//1)
    P2 = polarization_sum_massive(ni, np, k2, 1//1)
    e = tr * vtx * vtx_c * P1 * P2
    for _ in 1:10  # multi-pass until all factors are SP or Eps(all-Momentum)
        e = expand_scalar_product(contract(e))
        all(all(f isa FeynfeldX.Pair{Momentum,Momentum} || (f isa Eps && FeynfeldX._eps_all_momentum(f))
                for f in fk.factors) for (fk,_) in e.terms) && break
    end
    e
end

# ── Helper: fully contract fermion exchange (m_ν = 0) ──
function fermion_contract(chain_mom::DiracChain)
    tr = FeynfeldX._single_line_trace(DiracExpr(chain_mom))
    mi = LorentzIndex(:mu_k1, DimD()); mp = LorentzIndex(:mu_k1_, DimD())
    ni = LorentzIndex(:mu_k2, DimD()); np = LorentzIndex(:mu_k2_, DimD())
    P1 = polarization_sum_massive(mi, mp, k1, 1//1)
    P2 = polarization_sum_massive(ni, np, k2, 1//1)
    e = tr * P1 * P2
    for _ in 1:10
        e = expand_scalar_product(contract(e))
        all(all(f isa FeynfeldX.Pair{Momentum,Momentum} for f in fk.factors)
            for (fk,_) in e.terms) && break
    end
    e
end

# ── Build all symbolic expressions once ──
chain_sg, vtx_sg = amp_sg
chain_sZ, vtx_sZ = amp_sZ
chain_t_mom, _, _ = amp_t

sym_sg = gauge_contract(chain_sg, vtx_sg)
sym_sZ = gauge_contract(chain_sZ, vtx_sZ)
sym_t  = fermion_contract(chain_t_mom)
# Combined s-channel (s-γ + s-Z) — produces all diagonal + cross terms
combined_s = DiracExpr(vcat(chain_sg.terms, chain_sZ.terms))
sym_s_comb = gauge_contract(combined_s, vtx_sg)

# ── Stage A: diagonal terms at integer point ──
@testset "Stage A: diagonal |M_i|²" begin
    val_sg = _scalar(evaluate_sp(sym_sg; ctx=ctx_int))
    val_sZ = _scalar(evaluate_sp(sym_sZ; ctx=ctx_int))
    val_t  = _scalar(evaluate_sp(sym_t; ctx=ctx_int))
    @test val_sg isa Rational{Int} && val_sg != 0
    @test val_sZ isa Rational{Int} && val_sZ != 0
    @test val_t isa Rational{Int} && val_t != 0
    println("  T_sγ=$val_sg  T_sZ=$val_sZ  T_tν=$val_t")
end

# ── Stage A: diagonal-only integration ──
@testset "Stage A: diagonal sum vs Grozin" begin
    println("\n  √s/GeV   σ_diag [pb]   σ_Grozin [pb]   ratio")
    for sqrts in [170.0, 200.0, 300.0, 500.0]
        s_mw = sqrts^2 / const_MW^2; β = sqrt(1.0 - 4.0/s_mw)
        σ, _ = quadgk(-1.0, 1.0) do cosθ
            t_mw = -(s_mw/2)*(1-β*cosθ) + 1; u_mw = -(s_mw/2)*(1+β*cosθ) + 1
            sp = sp_values_2to2(Mandelstam(s_mw, t_mw, u_mw); m3_sq=1.0, m4_sq=1.0)
            m2 = 0.25 * (const_c_sg^2/s_mw^2 * evaluate_numeric(sym_sg, sp) +
                         const_c_sZ^2/(s_mw-const_MZ2)^2 * evaluate_numeric(sym_sZ, sp) +
                         const_c_t^2/t_mw^2 * evaluate_numeric(sym_t, sp))
            β/(32π*s_mw) * m2
        end
        σ_pb = σ/const_MW^2*const_PB; σ_G = sigma_ee_ww(sqrts^2)
        println("  $(lpad(Int(sqrts),4))    $(lpad(round(σ_pb,digits=3),9))    $(lpad(round(σ_G,digits=3),9))     $(round(σ_pb/σ_G,digits=3))")
        @test isfinite(σ_pb) && σ_pb > 0
    end
end

# ── Stage B: s-channel with γ-Z interference ──
# Ref: Denner1993, Eq. (11.9): M_s = M_sγ + M_sZ (coherent sum)
@testset "Stage B: s-channel γ-Z interference" begin
    # Cross-term = combined - diagonal_γ - diagonal_Z
    # _single_line_trace sums ALL (i,j) pairs, so cross includes both orientations.
    val_comb  = _scalar(evaluate_sp(sym_s_comb; ctx=ctx_int))
    val_sg_i  = _scalar(evaluate_sp(sym_sg; ctx=ctx_int))
    val_sZ_i  = _scalar(evaluate_sp(sym_sZ; ctx=ctx_int))
    val_cross = val_comb - val_sg_i - val_sZ_i
    @test val_cross != 0
    println("  T_cross(sγ×sZ) = $val_cross")

    println("\n  √s/GeV   σ(s+int) [pb]   σ_Grozin [pb]   ratio")
    for sqrts in [170.0, 200.0, 300.0, 500.0]
        s_mw = sqrts^2 / const_MW^2; β = sqrt(1.0 - 4.0/s_mw)
        σ, _ = quadgk(-1.0, 1.0) do cosθ
            t_mw = -(s_mw/2)*(1-β*cosθ) + 1; u_mw = -(s_mw/2)*(1+β*cosθ) + 1
            sp = sp_values_2to2(Mandelstam(s_mw, t_mw, u_mw); m3_sq=1.0, m4_sq=1.0)
            m2_sg = evaluate_numeric(sym_sg, sp)
            m2_sZ = evaluate_numeric(sym_sZ, sp)
            m2_comb = evaluate_numeric(sym_s_comb, sp)
            m2_cross = m2_comb - m2_sg - m2_sZ
            # w_γ² T_sγ + w_Z² T_sZ + w_γ w_Z T_cross (no factor of 2: T_cross has both orientations)
            w_γ = const_c_sg / s_mw; w_Z = const_c_sZ / (s_mw - const_MZ2)
            m2_s = w_γ^2*m2_sg + w_Z^2*m2_sZ + w_γ*w_Z*m2_cross
            m2 = 0.25 * (m2_s + const_c_t^2/t_mw^2 * evaluate_numeric(sym_t, sp))
            β/(32π*s_mw) * m2
        end
        σ_pb = σ/const_MW^2*const_PB; σ_G = sigma_ee_ww(sqrts^2)
        println("  $(lpad(Int(sqrts),4))     $(lpad(round(σ_pb,digits=3),9))    $(lpad(round(σ_G,digits=3),9))     $(round(σ_pb/σ_G,digits=3))")
        @test isfinite(σ_pb) && σ_pb > 0
    end
end

end  # testset
