# Layer 3b: Channel enumeration for 2->2 tree-level scattering.
#
# Key insight: tree-level 2->2 has exactly 3 topologies (s, t, u).
# Diagram generation is FILTERING (check vertex existence), not graph enumeration.
# No external library needed.
#
# Ref: SPIRAL_9_PLAN.md — "TreeDiagram22 is fully specified by
#      4 legs + 1 internal line + 2 vertices"

"""
    TreeChannel

A 2->2 tree-level scattering channel, specified by:
- `channel`: topology (:s, :t, :u)
- `exchanged`: field name of the internal propagator
- `legs`: (in1, in2, out1, out2) in standardized order
"""
struct TreeChannel
    channel::Symbol
    exchanged::Symbol
    legs::NTuple{4, ExternalLeg}
end

Base.show(io::IO, ch::TreeChannel) = print(io, "$(ch.channel)-channel [$(ch.exchanged)]")

"""
    vertex_legs(ch::TreeChannel)

Which external legs connect at each vertex, determined by channel topology.
Returns `((legL1, legL2), (legR1, legR2))`.
"""
function vertex_legs(ch::TreeChannel)
    _vertex_legs_by_type(ch.channel, ch.legs)
end

function _vertex_legs_by_type(ch_type::Symbol, legs::NTuple{4,ExternalLeg})
    in1, in2, out1, out2 = legs
    if ch_type == :s
        (in1, in2), (out1, out2)
    elseif ch_type == :t
        (in1, out1), (in2, out2)
    else  # :u
        (in1, out2), (in2, out1)
    end
end

"""
    tree_channels(model, rules, incoming, outgoing)

Enumerate all valid 2->2 tree-level channels by checking vertex existence
in the Feynman rules. Returns `Vector{TreeChannel}`.

Complexity: O(3 channels x N_fields).
"""
function tree_channels(model::AbstractModel, rules::FeynmanRules,
                       incoming::Vector{ExternalLeg}, outgoing::Vector{ExternalLeg})
    length(incoming) == 2 || error("tree_channels: need exactly 2 incoming")
    length(outgoing) == 2 || error("tree_channels: need exactly 2 outgoing")
    in1, in2 = incoming
    out1, out2 = outgoing
    legs = (in1, in2, out1, out2)
    channels = TreeChannel[]
    for ch_type in (:s, :t, :u)
        (l1, l2), (r1, r2) = _vertex_legs_by_type(ch_type, legs)
        for field in model_fields(model)
            _has_vertex(rules, l1.field_name, l2.field_name, field.name) || continue
            _has_vertex(rules, r1.field_name, r2.field_name, field.name) || continue
            # Boson exchange: fermion flow must be valid at fermion-pair vertices
            if field isa Field{Boson}
                if _both_fermion_legs(model, l1, l2)
                    _compatible_fermion_pair(l1, l2) || continue
                end
                if _both_fermion_legs(model, r1, r2)
                    _compatible_fermion_pair(r1, r2) || continue
                end
            end
            push!(channels, TreeChannel(ch_type, field.name, legs))
        end
    end
    channels
end

# Check if 3 fields form a valid vertex in any permutation.
function _has_vertex(rules::FeynmanRules, a::Symbol, b::Symbol, c::Symbol)
    for perm in ((a,b,c), (a,c,b), (b,a,c), (b,c,a), (c,a,b), (c,b,a))
        haskey(rules.vertices, perm) && return true
    end
    false
end

# Check if both external legs at a vertex are fermions.
function _both_fermion_legs(model::AbstractModel, a::ExternalLeg, b::ExternalLeg)
    get_field(model, a.field_name) isa Field{Fermion} &&
    get_field(model, b.field_name) isa Field{Fermion}
end

# Two fermion legs at a vertex must have compatible spinor positions:
# one bar (left) and one plain (right). Bar = outgoing particle or incoming antiparticle.
function _compatible_fermion_pair(a::ExternalLeg, b::ExternalLeg)
    bar_a = (a.incoming == a.antiparticle)  # incoming anti OR outgoing particle
    bar_b = (b.incoming == b.antiparticle)
    bar_a != bar_b  # one left, one right
end
