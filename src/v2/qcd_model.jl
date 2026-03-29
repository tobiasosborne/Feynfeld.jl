# QCD model: quark-gluon vertex (qqg) + triple gluon vertex (ggg).
#
# The qqg vertex has the same Lorentz structure as eeγ: γ^μ.
# The ggg vertex is momentum-dependent: V_{μ₁μ₂μ₃}(p1,p2,p3).
#
# Ref: Peskin & Schroeder, Eq. (16.10) for triple gauge vertex

struct QCDModel <: AbstractModel
    quark::Field{Fermion}
    gluon::Field{Boson}
end

model_name(::QCDModel) = :QCD
model_fields(m::QCDModel) = Field[m.quark, m.gluon]
gauge_groups(::QCDModel) = GaugeGroup[SU{3}()]
Base.show(io::IO, ::QCDModel) = print(io, "QCDModel(q, g)")

function qcd_model(; m_q=:zero)
    QCDModel(
        fermion(:q, m_q; charge=0//1),
        vector_boson(:g, :zero; self_conj=true),
    )
end

# Feynman rules: qqg + ggg vertices
function feynman_rules(model::QCDModel)
    gauge = first(gauge_groups(model))
    vertices = Dict{NTuple{3,Symbol}, VertexRule}()
    q, g = model.quark, model.gluon
    vertices[(q.name, q.name, g.name)] = VertexRule((q.name, q.name, g.name), :g_s)
    vertices[(g.name, g.name, g.name)] = VertexRule((g.name, g.name, g.name), :g_s)
    FeynmanRules(model, vertices, gauge)
end

# qqg vertex structure: γ^μ (same as QED)
function vertex_structure(::SU{N}, ::Fermion, ::Boson, mu::LorentzIndex) where N
    DiracExpr(DiracChain([DiracGamma(LISlot(mu))]))
end

"""
    triple_gauge_vertex(mu1, mu2, mu3, p1, p2, p3)

Triple gauge boson vertex V_{μ₁μ₂μ₃}(p1,p2,p3) with all momenta outgoing.

Ref: Peskin & Schroeder, Eq. (16.10)
"V^{abc}_{μνρ} = gf^{abc}[(k1-k2)_ρ g_{μν} + (k2-k3)_μ g_{νρ} + (k3-k1)_ν g_{μρ}]"
Returns AlgSum (coupling stripped).
"""
function triple_gauge_vertex(mu1::LorentzIndex, mu2::LorentzIndex, mu3::LorentzIndex,
                              p1::Momentum, p2::Momentum, p3::Momentum)
    _vtx_term(mu1, mu2, p1, p2, mu3) +
    _vtx_term(mu2, mu3, p2, p3, mu1) +
    _vtx_term(mu3, mu1, p3, p1, mu2)
end

# g_{ga,gb} × (pa - pb)_{gc} = Σ_m c_m pair(ga,gb) × pair(gc,m)
function _vtx_term(ga::LorentzIndex, gb::LorentzIndex,
                   pa::Momentum, pb::Momentum, gc::LorentzIndex)
    alg(pair(ga, gb)) * alg(pair(gc, pa)) - alg(pair(ga, gb)) * alg(pair(gc, pb))
end
