# Triple and quadruple gauge-boson vertices — pure Lorentz structure
# (model-agnostic).
#
# Lives in Layer 4 (not the QCD model) because these are plain algebraic
# helpers over AlgSum/Pair: used by the legacy gauge_exchange.jl path
# AND by the qgraf port's build_vertices (Phase 18b-4 ggg / Phase 18b-5
# gggg).

"""
    triple_gauge_vertex(mu1, mu2, mu3, p1, p2, p3)

Triple gauge boson vertex V_{μ₁μ₂μ₃}(p1,p2,p3) with all momenta outgoing.

Ref: refs/papers/PeskinSchroeder1995.djvu, Eq. (16.10)
"V^{abc}_{μνρ} = g f^{abc}[(k1-k2)_ρ g_{μν} + (k2-k3)_μ g_{νρ} + (k3-k1)_ν g_{μρ}]"

Returns an `AlgSum` (coupling g and colour factor f^{abc} stripped).
"""
const _MomLike = Union{Momentum, MomentumSum}

function triple_gauge_vertex(mu1::LorentzIndex, mu2::LorentzIndex, mu3::LorentzIndex,
                              p1::_MomLike, p2::_MomLike, p3::_MomLike)
    _vtx_term(mu1, mu2, p1, p2, mu3) +
    _vtx_term(mu2, mu3, p2, p3, mu1) +
    _vtx_term(mu3, mu1, p3, p1, mu2)
end

# g_{ga,gb} × (pa - pb)_{gc} — expands MomentumSum into individual pairs
function _vtx_term(ga::LorentzIndex, gb::LorentzIndex,
                   pa::_MomLike, pb::_MomLike, gc::LorentzIndex)
    g_ab = alg(pair(ga, gb))
    _mom_pair(gc, pa) * g_ab - _mom_pair(gc, pb) * g_ab
end

_mom_pair(li::LorentzIndex, p::Momentum) = alg(pair(li, p))
function _mom_pair(li::LorentzIndex, ms::MomentumSum)
    result = AlgSum()
    for (c, m) in ms.terms
        add!(result, alg(pair(li, m)), c)
    end
    result
end

"""
    quadruple_gauge_vertex(mu1, mu2, mu3, mu4)

Quadruple gauge boson contact vertex Lorentz tensor V_{μ₁μ₂μ₃μ₄}.

Ref: refs/papers/PeskinSchroeder1995.djvu, Eq. (16.5)-(16.6); cross-
checked against refs/FeynCalc/.../Feynman/GluonVertex.m:108-118 verbatim:
"gl4v = -I coup^2 ( SUNF[a,b,e] SUNF[c,d,e] (Pair[mu,la] Pair[nu,si]
                                            - Pair[mu,si] Pair[nu,la])
                 + SUNF[a,c,e] SUNF[b,d,e] (Pair[mu,nu] Pair[la,si]
                                            - Pair[mu,si] Pair[nu,la])
                 + SUNF[a,d,e] SUNF[b,c,e] (Pair[mu,nu] Pair[la,si]
                                            - Pair[mu,la] Pair[nu,si]))".

Returns the SUM of the three Lorentz tensors with each colour pair
factor `f^{xxe} f^{yye}` and the global `-i g_s²` prefactor stripped.
Colour algebra is bead `feynfeld-yewo`: until then, `AmplitudeBundle`
carries no colour field, and this helper produces the bare Lorentz
tensor that yewo will multiply by the per-pairing colour structure.
"""
function quadruple_gauge_vertex(mu1::LorentzIndex, mu2::LorentzIndex,
                                 mu3::LorentzIndex, mu4::LorentzIndex)
    g(a, b) = alg(pair(a, b))
    # T_{(12)(34)} : colour pairing (a,b|c,d) → (g_{μ₁μ₃}g_{μ₂μ₄} − g_{μ₁μ₄}g_{μ₂μ₃})
    T12_34 = g(mu1, mu3) * g(mu2, mu4) - g(mu1, mu4) * g(mu2, mu3)
    # T_{(13)(24)} : colour pairing (a,c|b,d) → (g_{μ₁μ₂}g_{μ₃μ₄} − g_{μ₁μ₄}g_{μ₂μ₃})
    T13_24 = g(mu1, mu2) * g(mu3, mu4) - g(mu1, mu4) * g(mu2, mu3)
    # T_{(14)(23)} : colour pairing (a,d|b,c) → (g_{μ₁μ₂}g_{μ₃μ₄} − g_{μ₁μ₃}g_{μ₂μ₄})
    T14_23 = g(mu1, mu2) * g(mu3, mu4) - g(mu1, mu3) * g(mu2, mu4)
    T12_34 + T13_24 + T14_23
end
