# Triple gauge-boson vertex — pure Lorentz structure (model-agnostic).
#
# Lives in Layer 4 (not the QCD model) because it is a plain algebraic
# helper over AlgSum/Pair: it is used by the legacy gauge_exchange.jl
# path AND by the qgraf port's build_vertices (Phase 18b-4, ggg vertex).

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
