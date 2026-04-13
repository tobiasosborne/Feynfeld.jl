# Electroweak Standard Model: field table + vertex dispatch.
#
# NOT a Lagrangian parser — explicit vertex table for the SM interactions
# needed by tree-level e+e- -> W+W- (3 diagrams: s-gamma, s-Z, t-nu_e).
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.1)-(11.2)
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Table 10.3

struct EWModel <: AbstractModel
    electron::Field{Fermion}
    neutrino::Field{Fermion}
    photon::Field{Boson}
    z_boson::Field{Boson}
    w_boson::Field{Boson}
end

model_name(::EWModel) = :EW
function model_fields(m::EWModel)
    Field[m.electron, m.neutrino, m.photon, m.z_boson, m.w_boson]
end
gauge_groups(::EWModel) = GaugeGroup[SU{2}(), U1()]
Base.show(io::IO, ::EWModel) = print(io, "EWModel(e, ν_e, γ, Z, W)")

function ew_model()
    EWModel(
        fermion(:e, :zero; charge=-1//1),
        fermion(:nu_e, :zero; charge=0//1),
        vector_boson(:gamma, :zero),
        vector_boson(:Z, :M_Z),
        vector_boson(:W, :M_W; self_conj=false),
    )
end

# Feynman rules: 5 vertex types
function feynman_rules(model::EWModel)
    gauge = U1()  # gauge parameter only used for vertex_structure dispatch
    vertices = Dict{Tuple, VertexRule}()
    # eeγ: standard QED
    vertices[(:e, :e, :gamma)] = VertexRule((:e, :e, :gamma), :e)
    # eeZ: chiral coupling (g_V - g_A γ5)γ^μ
    vertices[(:e, :e, :Z)] = VertexRule((:e, :e, :Z), :e_Z)
    # eνW: charged current (1-γ5)/2 γ^μ
    vertices[(:e, :nu_e, :W)] = VertexRule((:e, :nu_e, :W), :g_W)
    # WWγ: triple gauge
    vertices[(:W, :W, :gamma)] = VertexRule((:W, :W, :gamma), :e)
    # WWZ: triple gauge
    vertices[(:W, :W, :Z)] = VertexRule((:W, :W, :Z), :g_WWZ)
    FeynmanRules(model, vertices, gauge)
end
