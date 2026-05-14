# QCD model: quark-gluon vertex (qqg) + triple gluon vertex (ggg).
#
# The qqg vertex has the same Lorentz structure as eeγ: γ^μ.
# The ggg vertex is momentum-dependent: V_{μ₁μ₂μ₃}(p1,p2,p3).
#
# Ref: Peskin & Schroeder, Eq. (16.10) for triple gauge vertex

#  Ghost is anticommuting (fermion statistics) but Lorentz-scalar; modelled
#  as Field{Fermion} for diagram counting. Layer-4 Lorentz refinement TBD.
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd "[ghost,antighost,-1]"
struct QCDModel <: AbstractModel
    quark::Field{Fermion}
    gluon::Field{Boson}
    ghost::Field{Fermion}
end

model_name(::QCDModel) = :QCD
model_fields(m::QCDModel) = Field[m.quark, m.gluon, m.ghost]
gauge_groups(::QCDModel) = GaugeGroup[SU{3}()]
Base.show(io::IO, ::QCDModel) = print(io, "QCDModel(q, g, ghost)")

function qcd_model(; m_q=:zero)
    QCDModel(
        fermion(:q, m_q; charge=0//1),
        vector_boson(:g, :zero; self_conj=true),
        fermion(:ghost, :zero; charge=0//1),
    )
end

# Feynman rules: qqg + ggg vertices
function feynman_rules(model::QCDModel)
    gauge = first(gauge_groups(model))
    vertices = Dict{Tuple, VertexRule}()
    q, g = model.quark, model.gluon
    vertices[(q.name, q.name, g.name)] = VertexRule((q.name, q.name, g.name), :g_s)
    vertices[(g.name, g.name, g.name)] = VertexRule((g.name, g.name, g.name), :g_s)
    # 4-gluon contact vertex.
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd "[gluon,gluon,gluon,gluon]"
    vertices[(g.name, g.name, g.name, g.name)] =
        VertexRule((g.name, g.name, g.name, g.name), :g_s)
    # Ghost-gluon vertex.  Vertex tuple uses bare names; _expand_vertex maps
    # the fermion pair (ghost, ghost) → (ghost, ghost_bar) downstream.
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd "[antighost,ghost,gluon]"
    gh = model.ghost
    vertices[(gh.name, gh.name, g.name)] = VertexRule((gh.name, gh.name, g.name), :g_s)
    FeynmanRules(model, vertices, gauge)
end

# `triple_gauge_vertex` (the ggg / WWγ Lorentz structure) lives in
# src/v2/gauge_vertex.jl — it is a model-agnostic Layer-4 algebraic
# helper shared by gauge_exchange.jl and the qgraf port's build_vertices.
