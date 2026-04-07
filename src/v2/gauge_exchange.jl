# Gauge boson exchange amplitude: fermion pair ↔ triple gauge vertex.
#
# qq̄→gg s-channel: one vertex is qqg (fermion line), the other is ggg.
# Returns (chain::DiracExpr, vertex::AlgSum):
#   chain = ū(p) γ^ρ v(p')  (fermion line with propagator index ρ)
#   vertex = V_{ρ,μ₁,μ₂}(k₁, k₂, -q)  (triple gauge, all-outgoing convention)
#
# Ref: Peskin & Schroeder, Eq. (16.10) for triple gauge vertex

function _build_gauge_exchange(ch::TreeChannel, rules::FeynmanRules, model::AbstractModel)
    (legL1, legL2), (legR1, legR2) = vertex_legs(ch)

    # Identify fermion pair vs boson pair
    if _both_fermion_legs(model, legL1, legL2)
        ferm1, ferm2 = legL1, legL2
        bos1, bos2 = legR1, legR2
    else
        ferm1, ferm2 = legR1, legR2
        bos1, bos2 = legL1, legL2
    end

    # Propagator Lorentz index (shared between fermion vertex and gauge vertex)
    rho = LorentzIndex(Symbol(:rho_, ch.channel), DimD())

    # Fermion chain with vertex structure from rules
    vtx = _lookup_vertex(rules, ferm1, ferm2, ch.exchanged, rho)
    chain = _fermion_line_chain(ferm1, ferm2, vtx)

    # Triple gauge vertex in all-outgoing convention:
    # propagator carries q into the vertex → outgoing momentum is -q
    q = propagator_momentum(ch)
    neg_q = MomentumSum([(-c, m) for (c, m) in q.terms])
    mu1 = LorentzIndex(Symbol(:mu_, bos1.momentum.name), DimD())
    mu2 = LorentzIndex(Symbol(:mu_, bos2.momentum.name), DimD())
    vtx = triple_gauge_vertex(rho, mu1, mu2, neg_q, bos1.momentum, bos2.momentum)

    (chain, vtx)
end
