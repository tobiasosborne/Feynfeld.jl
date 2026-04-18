# ee→WW multi-channel validation against Grozin formula.
#
# Ground truth:
#   Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.9)-(11.10)
#   Ref: refs/FeynCalc/.../AnelEl-WW.m, lines 80, 95-98, 193-198
#   Ref: src/v2/ew_cross_section.jl (Grozin formula, verified vs Denner Tab. 11.4)
#
# Coupling constants NOT in pipeline output (must multiply explicitly):
#   Ref: Denner1993, Eq. (11.9); PDG2024, Table 10.3
#   c_sγ = +e² = +4πα                (eeγ × WWγ, C_γWW = +1)
#   c_sZ = -e²/(2sin²θ_W)           (eeZ × WWZ, C_ZWW = -cW/sW → sign from Denner A.7)
#   c_tν = +e²/(2sin²θ_W)           (eνW × eνW, P_L already in vertex_structure)
# The SIGN of c_sZ comes from the VVV coupling C: γWW has C=+1, ZWW has C=-cW/sW.
# The magnitude |cW/sW| is absorbed into the coupling product e²/(2sin²θ_W).
# The relative phase matters for s×t and γ×Z interference terms.

using Test
using Feynfeld
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

# Coupling constants — Ref: Denner1993, Eqs. (11.9), (A.7)
const_e2   = 4π * EW_ALPHA
const_c_sg = const_e2                                             # γWW: C = +1
const_c_sZ = gauge_coupling_phase(Val(:g_WWZ)) * const_e2 / (2 * EW_SIN2_W)  # ZWW: C = -cW/sW
const_c_t  = const_e2 / (2 * EW_SIN2_W)                          # eνW × eνW
const_MZ2  = 1.0 / (1.0 - EW_SIN2_W)  # M_Z² in M_W units
const_MW   = EW_M_W
const_PB   = 3.894e8  # GeV⁻² → pb

# Shared integer kinematic point (M_W units)
s_val = 8//1; t_val = -1//1; u_val = 4//1 - s_val - t_val
ctx_int = sp_context(
    (:p1,:p1)=>0//1, (:p2,:p2)=>0//1, (:k1,:k1)=>1//1, (:k2,:k2)=>1//1,
    (:p1,:p2)=>s_val//2, (:p1,:k1)=>-t_val//2, (:p1,:k2)=>-u_val//2,
    (:p2,:k1)=>-u_val//2, (:p2,:k2)=>-t_val//2, (:k1,:k2)=>(s_val-2)//2)

# Helper: extract scalar coefficient from evaluated AlgSum.
# Only the empty-FactorKey term is the true scalar; Eps terms are parity-odd artifacts.
function _scalar(r)
    scalar_key = FactorKey()
    haskey(r.terms, scalar_key) || return 0//1
    c = r.terms[scalar_key]
    c isa DimPoly ? evaluate_dim(c) : c
end

# ── Helper: fully contract gauge exchange expression ──
function gauge_contract(chain_de::DiracExpr, vtx::AlgSum)
    tr = Feynfeld._single_line_trace(chain_de)
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
        all(all(f isa Feynfeld.Pair{Momentum,Momentum} || (f isa Eps && Feynfeld._eps_all_momentum(f))
                for f in fk.factors) for (fk,_) in e.terms) && break
    end
    e
end

# ── Helper: fully contract fermion exchange (m_ν = 0) ──
function fermion_contract(chain_mom::DiracChain)
    tr = Feynfeld._single_line_trace(DiracExpr(chain_mom))
    mi = LorentzIndex(:mu_k1, DimD()); mp = LorentzIndex(:mu_k1_, DimD())
    ni = LorentzIndex(:mu_k2, DimD()); np = LorentzIndex(:mu_k2_, DimD())
    P1 = polarization_sum_massive(mi, mp, k1, 1//1)
    P2 = polarization_sum_massive(ni, np, k2, 1//1)
    e = tr * P1 * P2
    for _ in 1:10
        e = expand_scalar_product(contract(e))
        all(all(f isa Feynfeld.Pair{Momentum,Momentum} || (f isa Eps && Feynfeld._eps_all_momentum(f))
                for f in fk.factors) for (fk,_) in e.terms) && break
    end
    e
end

# ── Helper: cross-trace for one orientation (s-fwd, t-conj) ──
# Computes Σ_{i∈fwd, j∈conj} c_i c_j sign_j Tr[compl_R Γ̃_j compl_L Γ_i]
# Uses the same internal machinery as _single_line_trace (spin_sum.jl:88-98).
function cross_trace_one_way(de_fwd::DiracExpr, de_conj::DiracExpr)
    sp_L = first(de_fwd.terms[1][2].elements)::Spinor
    sp_R = last(de_fwd.terms[1][2].elements)::Spinor
    _, mass_R, mom_R = Feynfeld._completeness(sp_R)
    _, mass_L, mom_L = Feynfeld._completeness(sp_L)
    result = AlgSum()
    for (ci, chain_i) in de_fwd.terms
        gammas_i = DiracGamma[e for e in chain_i.elements[2:end-1]]
        for (cj, chain_j) in de_conj.terms
            gammas_j = DiracGamma[e for e in chain_j.elements[2:end-1]]
            sign_j = Feynfeld._conj_gamma5_sign(gammas_j)
            gammas_conj_j = Feynfeld._conjugate_gammas(gammas_j)
            tr = Feynfeld._build_and_trace(mom_R, mass_R, gammas_conj_j,
                                            mom_L, mass_L, gammas_i)
            result = result + ci * cj * sign_j * tr
        end
    end
    result
end

# Multi-pass contraction (shared logic)
function full_contract(e::AlgSum)
    for _ in 1:10
        e = expand_scalar_product(contract(e))
        all(all(f isa Feynfeld.Pair{Momentum,Momentum} ||
               (f isa Eps && Feynfeld._eps_all_momentum(f))
            for f in fk.factors) for (fk,_) in e.terms) && break
    end
    e
end

# ── Helper: cross contraction for s × t interference ──
# Returns (fwd_expr, conj_expr) — the two orientations of the cross-term.
function cross_contract_st(de_s::DiracExpr, vtx_s::AlgSum, chain_t_mom::DiracChain)
    de_t = DiracExpr([(alg(1), chain_t_mom)])
    neg_q = MomentumSum([(-1//1, p1), (-1//1, p2)])
    mi = LorentzIndex(:mu_k1, DimD()); mp = LorentzIndex(:mu_k1_, DimD())
    ni = LorentzIndex(:mu_k2, DimD()); np = LorentzIndex(:mu_k2_, DimD())
    P1 = polarization_sum_massive(mi, mp, k1, 1//1)
    P2 = polarization_sum_massive(ni, np, k2, 1//1)

    # Orientation 1: s-fwd × t-conj
    # Free indices: rho_s (s-channel) + mu_k1_, mu_k2_ (conjugated t-channel)
    # Contract with vtx_s(rho_s, mu_k1, mu_k2) × P1(mu_k1, mu_k1_) × P2(mu_k2, mu_k2_)
    tr_fwd = cross_trace_one_way(de_s, de_t)
    e_fwd = full_contract(tr_fwd * vtx_s * P1 * P2)

    # Orientation 2: t-fwd × s-conj
    # Free indices: mu_k1, mu_k2 (t-channel) + rho_s_ (conjugated s-channel)
    # Contract with vtx_c(rho_s_, mu_k1_, mu_k2_) × P1 × P2
    vtx_c = triple_gauge_vertex(LorentzIndex(:rho_s_, DimD()),
        LorentzIndex(:mu_k1_, DimD()), LorentzIndex(:mu_k2_, DimD()), neg_q, k1, k2)
    tr_conj = cross_trace_one_way(de_t, de_s)
    e_conj = full_contract(tr_conj * vtx_c * P1 * P2)

    (e_fwd, e_conj)
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

# Stage C: s × t cross-terms (both orientations)
sym_sg_t_fwd, sym_sg_t_conj = cross_contract_st(chain_sg, vtx_sg, chain_t_mom)
sym_sZ_t_fwd, sym_sZ_t_conj = cross_contract_st(chain_sZ, vtx_sZ, chain_t_mom)

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

# ── Stage C: full cross-section with s×t gauge cancellation ──
# Ref: Denner1993, Eqs. (11.16)-(11.17): gauge cancellation makes σ ~ log(s)/s
# Without s×t terms, σ grows as s²/M_W⁴ (unphysical).
@testset "Stage C: full cross-section vs Grozin" begin
    # Sanity: cross-terms are non-zero at integer kinematics
    val_sg_t_f = _scalar(evaluate_sp(sym_sg_t_fwd; ctx=ctx_int))
    val_sg_t_c = _scalar(evaluate_sp(sym_sg_t_conj; ctx=ctx_int))
    @test val_sg_t_f != 0 || val_sg_t_c != 0

    println("\n  √s/GeV   σ_full [pb]   σ_Grozin [pb]   ratio")
    for sqrts in [170.0, 200.0, 300.0, 500.0]
        s_mw = sqrts^2 / const_MW^2; β = sqrt(1.0 - 4.0/s_mw)
        σ, _ = quadgk(-1.0, 1.0; rtol=1e-6) do cosθ
            t_mw = -(s_mw/2)*(1-β*cosθ) + 1; u_mw = -(s_mw/2)*(1+β*cosθ) + 1
            sp = sp_values_2to2(Mandelstam(s_mw, t_mw, u_mw); m3_sq=1.0, m4_sq=1.0)

            # Diagonal terms
            m2_sg = evaluate_numeric(sym_sg, sp)
            m2_sZ = evaluate_numeric(sym_sZ, sp)
            m2_t  = evaluate_numeric(sym_t, sp)

            # s-channel γ-Z interference (Stage B)
            m2_s_comb = evaluate_numeric(sym_s_comb, sp)
            m2_sgZ_cross = m2_s_comb - m2_sg - m2_sZ

            # s × t cross-terms (Stage C)
            m2_sg_t = evaluate_numeric(sym_sg_t_fwd, sp) + evaluate_numeric(sym_sg_t_conj, sp)
            m2_sZ_t = evaluate_numeric(sym_sZ_t_fwd, sp) + evaluate_numeric(sym_sZ_t_conj, sp)

            # Propagator weights
            w_γ = const_c_sg / s_mw
            w_Z = const_c_sZ / (s_mw - const_MZ2)
            w_t = const_c_t / t_mw

            # Full |M|² / 4 (spin average)
            m2 = 0.25 * (
                w_γ^2 * m2_sg + w_Z^2 * m2_sZ + w_γ * w_Z * m2_sgZ_cross +
                w_t^2 * m2_t +
                w_γ * w_t * m2_sg_t + w_Z * w_t * m2_sZ_t
            )
            β/(32π*s_mw) * m2
        end
        σ_pb = σ/const_MW^2*const_PB; σ_G = sigma_ee_ww(sqrts^2)
        ratio = σ_pb / σ_G
        println("  $(lpad(Int(sqrts),4))    $(lpad(round(σ_pb,digits=3),9))    $(lpad(round(σ_G,digits=3),9))     $(round(ratio,digits=4))")
        # Full gauge cancellation: exact match at all energies.
        # Ref: Denner1993 Eq. (11.17): σ ~ log(s)/s at high energy.
        tol = 1e-4
        @test abs(ratio - 1.0) < tol
    end
end

end  # testset
