# Tensor Integral Decomposition (TID): evaluate 1-loop box integrals.
#
# Given the traced AlgSum (with SP(q,p_i), SP(q,q), and external SP terms),
# numerically evaluate ∫ d^D q N(q)/[D₀D₁D₂D₃].
#
# Rank 0: coefficient × D₀
# Rank 1: ∫ q·p_a / [box] = Σ_i (K_i·p_a) D_i
# Rank 2: PV cancellation → rank-1 triangle + rank-1 box
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.7)-(4.18)

"""
    evaluate_box_integral(numerator, sp_vals, denoms, s, man_inv; mu2)

Evaluate the loop integral for a box diagram.

- `numerator`: AlgSum with SP(q,...) from interference trace
- `sp_vals`: external kinematics as Dict{(Symbol,Symbol), Float64}
- `denoms`: BoxDenominators with accumulated momenta
- `s`: CM energy squared
- `man_inv`: the Mandelstam invariant for p₁₃ (t for direct, u for crossed)
"""
function evaluate_box_integral(numerator::AlgSum,
                               sp_vals::Dict{Tuple{Symbol,Symbol}, Float64},
                               denoms::BoxDenominators,
                               s::Float64, man_inv::Float64;
                               mu2::Float64 = 1.0)
    invs = (0.0, 0.0, 0.0, 0.0, s, man_inv)
    masses = (0.0, 0.0, 0.0, 0.0)
    K1, K2, K3 = denoms.accum_mom_names

    # Pre-compute scalar and tensor D functions
    d0 = _D0_evaluate(invs..., masses...)
    d1 = _eval_D_tensor(D1(invs..., masses...); mu2=mu2)
    d2 = _eval_D_tensor(D2(invs..., masses...); mu2=mu2)
    d3 = _eval_D_tensor(D3(invs..., masses...); mu2=mu2)

    # Pre-compute sub-triangle C₀ for SP(q,q) rank-2
    # tri_0 (remove prop 0): C₀(p12,p23,p13, m1²,m2²,m3²)
    c0_tri0 = _C0_analytical(0.0, 0.0, man_inv, 0.0, 0.0, 0.0)

    # f_j = m_j² - m₀² - K_j² (massless: f_j = -K_j²)
    f = (0.0, -s, 0.0)  # f₁=-K₁²=0, f₂=-K₂²=-s, f₃=-K₃²=0

    total = ComplexF64(0)
    for (fk, coeff) in numerator.terms
        c_val = Float64(evaluate_dim(coeff))
        iszero(c_val) && continue

        q_sps = ScalarProduct[]
        ext_val = 1.0
        for fac in fk.factors
            if _involves_q(fac)
                push!(q_sps, fac)
            else
                ext_val *= _eval_factor(fac, sp_vals)
            end
        end
        iszero(ext_val) && continue
        c_ext = c_val * ext_val

        # SP(q,q) counts as rank 2 (q²), not rank 1
        rank = sum(sp -> (sp.a.name == :q && sp.b.name == :q) ? 2 : 1, q_sps; init=0)
        if rank == 0
            total += c_ext * d0
        elseif rank == 1
            pa = _q_partner(q_sps[1])
            total += c_ext * _rank1(pa, K1, K2, K3, d1, d2, d3, sp_vals)
        elseif rank == 2
            total += c_ext * _rank2(q_sps, sp_vals, K1, K2, K3,
                                     d1, d2, d3, c0_tri0, f, invs, masses; mu2=mu2)
        else
            error("TID rank $rank > 2")
        end
    end
    total
end

# ---- Helpers ----
_involves_q(f::ScalarProduct) = (f.a.name == :q || f.b.name == :q)
_involves_q(::Any) = false

function _q_partner(sp::ScalarProduct)
    sp.a.name == :q && return sp.b.name
    sp.b.name == :q && return sp.a.name
    error("SP does not involve q")
end

# K · p_a where K is Momentum or MomentumSum
_Kdot(K::Momentum, pa::Symbol, sv) = _splookup(sv, K.name, pa)
function _Kdot(K::MomentumSum, pa::Symbol, sv)
    sum(Float64(c) * _splookup(sv, m.name, pa) for (c, m) in K.terms)
end

function _splookup(sv::Dict{Tuple{Symbol,Symbol},Float64}, a::Symbol, b::Symbol)
    key = a <= b ? (a, b) : (b, a)
    haskey(sv, key) && return sv[key]
    haskey(sv, (b, a)) && return sv[(b, a)]
    error("SP($a,$b) not in sp_vals")
end

# ---- Rank 1 ----
# ∫ q·p_a / [D₀...D₃] = Σ_i (K_i · p_a) D_i
# Ref: Denner1993 Eq. (4.13)
function _rank1(pa, K1, K2, K3, d1, d2, d3, sv)
    _Kdot(K1, pa, sv) * d1 + _Kdot(K2, pa, sv) * d2 + _Kdot(K3, pa, sv) * d3
end

# ---- Rank 2 ----
# SP(q,q): ∫ q²/[box] = C₀⁽⁰⁾ + m₀²D₀  (massless: just C₀⁽⁰⁾)
# SP(q,pa)×SP(q,pb): PV cancellation → rank-1 triangle + box integrals.
#   Uses 2q·K_j = D_j - D₀ + f_j to cancel propagators.
#   Sub-triangle C₁/C₂ now computable: B₀(0,0,0) handled by COLLIER.
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (4.7)-(4.9)
function _rank2(q_sps, sv, K1, K2, K3,
                d1, d2, d3, c0_tri0, f, invs, masses; mu2=1.0)
    p10, p12, p23, p30, p20, p13 = invs

    # SP(q,q): single FactorKey with both slots = :q
    if length(q_sps) == 1
        @assert q_sps[1].a.name == :q && q_sps[1].b.name == :q
        return c0_tri0 + masses[1] * _D0_evaluate(invs..., masses...)
    end

    # SP(q,pa) × SP(q,pb): decompose pa in K basis, then PV cancel
    pa = _q_partner(q_sps[1])
    pb = _q_partner(q_sps[2])

    # Decompose pa = α₁K₁ + α₂K₂ + α₃K₃ via Gram matrix
    Ks = (K1, K2, K3)
    b = ntuple(j -> _Kdot(Ks[j], pa, sv), 3)
    g = ntuple(6) do idx  # g11,g12,g13,g22,g23,g33
        ((i,j) = ((1,1),(1,2),(1,3),(2,2),(2,3),(3,3))[idx]; _Kdot(Ks[i], Ks[j], sv))
    end
    g11,g12,g13,g22,g23,g33 = g
    det_G = g11*(g22*g33-g23^2) - g12*(g12*g33-g23*g13) + g13*(g12*g23-g22*g13)
    abs(det_G) < 1e-30 && error("TID rank-2: Gram singular (det=$det_G)")
    alphas = (
        ((g22*g33-g23^2)*b[1] + (g13*g23-g12*g33)*b[2] + (g12*g23-g13*g22)*b[3]) / det_G,
        ((g23*g13-g12*g33)*b[1] + (g11*g33-g13^2)*b[2] + (g12*g13-g11*g23)*b[3]) / det_G,
        ((g12*g23-g13*g22)*b[1] + (g13*g12-g11*g23)*b[2] + (g11*g22-g12^2)*b[3]) / det_G,
    )

    # Sub-triangle invariants (Denner d_tensor.jl pattern):
    # tri_j: remove prop j from box → triangle of remaining 3 props
    tri_invs = ((p12,p23,p13), (p20,p23,p30), (p10,p13,p30), (p10,p12,p20))
    tri_masses_arr = (
        (masses[2],masses[3],masses[4]),  # tri_0: props 1,2,3
        (masses[1],masses[3],masses[4]),  # tri_1: props 0,2,3
        (masses[1],masses[2],masses[4]),  # tri_2: props 0,1,3
        (masses[1],masses[2],masses[3]),  # tri_3: props 0,1,2
    )
    # Sub-triangle accumulated momenta K'₁, K'₂:
    # tri_0: K'₁=K₂-K₁, K'₂=K₃-K₁
    # tri_j(j≥1): K'₁,K'₂ are subsets of {K₁,K₂,K₃} minus K_j index

    box_r1_pb = _rank1(pb, K1, K2, K3, d1, d2, d3, sv)
    result = ComplexF64(0)
    for (j, aj) in enumerate(alphas)
        abs(aj) < 1e-15 && continue
        tri_j = _tri_rank1(j, pb, sv, K1, K2, K3, tri_invs[j+1], tri_masses_arr[j+1]; mu2)
        tri_0 = _tri_rank1(0, pb, sv, K1, K2, K3, tri_invs[1], tri_masses_arr[1]; mu2)
        result += aj * 0.5 * (tri_j - tri_0 + f[j] * box_r1_pb)
    end
    result
end

# Rank-1 sub-triangle integral: ∫ q·pb / [tri_j]
# = (K'₁·pb) C₁ + (K'₂·pb) C₂
# Ref: Denner1993 Eq. (4.6)-(4.8) for C-tensor reduction
function _tri_rank1(j::Int, pb::Symbol, sv, K1, K2, K3, tri_inv, tri_mass; mu2=1.0)
    c1 = evaluate(C1(tri_inv..., tri_mass...); mu2=mu2)
    c2 = evaluate(C2(tri_inv..., tri_mass...); mu2=mu2)
    if j == 0  # K'₁=K₂-K₁, K'₂=K₃-K₁
        k1p = _Kdot(K2, pb, sv) - _Kdot(K1, pb, sv)
        k2p = _Kdot(K3, pb, sv) - _Kdot(K1, pb, sv)
    elseif j == 1  # K'₁=K₂, K'₂=K₃
        k1p = _Kdot(K2, pb, sv); k2p = _Kdot(K3, pb, sv)
    elseif j == 2  # K'₁=K₁, K'₂=K₃
        k1p = _Kdot(K1, pb, sv); k2p = _Kdot(K3, pb, sv)
    else  # j==3: K'₁=K₁, K'₂=K₂
        k1p = _Kdot(K1, pb, sv); k2p = _Kdot(K2, pb, sv)
    end
    k1p * c1 + k2p * c2
end

# Generalized _Kdot for K dotted with another Momentum/MomentumSum
_Kdot(K::Momentum, K_ref::Momentum, sv) = _splookup(sv, K.name, K_ref.name)
function _Kdot(K::Momentum, K_ref::MomentumSum, sv)
    sum(Float64(c) * _splookup(sv, K.name, m.name) for (c, m) in K_ref.terms)
end
function _Kdot(K::MomentumSum, K_ref::Momentum, sv)
    sum(Float64(c) * _splookup(sv, m.name, K_ref.name) for (c, m) in K.terms)
end
function _Kdot(K::MomentumSum, K_ref::MomentumSum, sv)
    val = 0.0
    for (ci, mi) in K.terms, (cj, mj) in K_ref.terms
        val += Float64(ci) * Float64(cj) * _splookup(sv, mi.name, mj.name)
    end
    val
end
