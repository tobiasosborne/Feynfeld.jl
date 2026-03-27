# Layer 3: Feynman diagrams — callable amplitude builders.
#
# Design patterns:
#   - Diagram as callable struct (diagram(momenta) → amplitude)
#   - Topology-first: define topology, then dress with fields
#   - Hard-coded tree topologies for tracer bullet; full generation deferred

struct ExternalLeg
    field_name::Symbol
    momentum::Momentum
    incoming::Bool
    antiparticle::Bool
end

struct FeynmanDiagram
    name::Symbol
    external::Vector{ExternalLeg}
    topology::Symbol  # :s_channel, :t_channel, :u_channel, :triangle, etc.
end

function Base.show(io::IO, d::FeynmanDiagram)
    in_str = join([l.field_name for l in d.external if l.incoming], ",")
    out_str = join([l.field_name for l in d.external if !l.incoming], ",")
    print(io, "$(d.name): $(in_str) → $(out_str) [$(d.topology)]")
end

# ---- Diagram is callable: build the amplitude expression ----
# For e+e- → μ+μ- s-channel, returns the two DiracChains
# that form the amplitude (one per fermion line).
function build_amplitude(d::FeynmanDiagram, rules::FeynmanRules)
    if d.topology == :s_channel_ee_mumu
        return _build_ee_mumu_s(d, rules)
    end
    error("Unknown topology: $(d.topology)")
end

function _build_ee_mumu_s(d::FeynmanDiagram, rules::FeynmanRules)
    # e+e- → γ* → μ+μ-
    # M = (ie)² / s × [v̄(p2) γ^μ u(p1)] × [ū(k1) γ_μ v(k2)]
    #
    # We return the two Dirac chains. The photon propagator
    # (-g_{μν}/s) contracts μ↔ν. The coupling e² and 1/s are
    # overall factors tracked separately.

    p1 = d.external[1].momentum  # e-
    p2 = d.external[2].momentum  # e+
    k1 = d.external[3].momentum  # μ-
    k2 = d.external[4].momentum  # μ+

    mu = LorentzIndex(:mu)

    chain_e = dot(vbar(p2), GA(:mu), u(p1))
    chain_mu = dot(ubar(k1), GA(:mu), v(k2))

    (chain_e, chain_mu)
end

# ---- Tree diagram generator for e+e- → μ+μ- ----
function tree_diagrams(model::AbstractModel, incoming, outgoing)
    mn = model_name(model)
    if mn == :QED && _is_ee_mumu(incoming, outgoing)
        return _tree_ee_mumu(incoming, outgoing)
    end
    error("Diagram generation for $incoming → $outgoing not implemented")
end

function _is_ee_mumu(incoming, outgoing)
    in_names = Set(l.field_name for l in incoming)
    out_names = Set(l.field_name for l in outgoing)
    in_names == Set([:e, :e]) && out_names == Set([:mu, :mu])
end

function _tree_ee_mumu(incoming, outgoing)
    [FeynmanDiagram(:ee_mumu_tree, [incoming; outgoing], :s_channel_ee_mumu)]
end
