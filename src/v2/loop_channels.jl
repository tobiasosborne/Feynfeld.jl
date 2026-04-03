# Layer 3b: 1-loop box channel enumeration for 2→2 scattering.
#
# Extends channels.jl (tree-level) to 1-loop box topologies.
# Key insight: a 2→2 QED box has exactly 2 topologies (direct + crossed).
# The loop has 4 internal propagators: 2 photon + 2 fermion.
#
# In a box, each fermion line connects a particle-antiparticle PAIR of the
# SAME flavor (e.g. e⁻-e⁺). The two photons connect the electron line to
# the muon line. "Direct" vs "crossed" refers to which photon connects
# which vertices on the two lines.
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 7
# Ref: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/ElAel-MuAmu.m

"""
    LoopChannel

A 1-loop box channel for 2→2 scattering, specified by:
- `topology`: `:direct_box` or `:crossed_box`
- `internal_fields`: field names of 4 internal propagators
- `legs`: (in1, in2, out1, out2) in standardized order
- `loop_momentum`: name of the loop integration variable
"""
struct LoopChannel
    topology::Symbol
    internal_fields::NTuple{4, Symbol}
    legs::NTuple{4, ExternalLeg}
    loop_momentum::Symbol
end

Base.show(io::IO, ch::LoopChannel) = print(io, "$(ch.topology) [$(join(ch.internal_fields, ","))]")

"""
    box_channels(model, rules, incoming, outgoing)

Enumerate all valid 1-loop box topologies for 2→2 scattering.
A box has two fermion lines (same-flavor pairs) connected by two virtual bosons.

For e⁻e⁺ → μ⁻μ⁺:
- Fermion line 1: (e⁻, e⁺) — electron line
- Fermion line 2: (μ⁻, μ⁺) — muon line
- Direct box: photon connects e⁻ vertex to μ⁻ vertex
- Crossed box: photon connects e⁻ vertex to μ⁺ vertex
"""
function box_channels(model::AbstractModel, rules::FeynmanRules,
                      incoming::Vector{ExternalLeg}, outgoing::Vector{ExternalLeg})
    length(incoming) == 2 || error("box_channels: need exactly 2 incoming")
    length(outgoing) == 2 || error("box_channels: need exactly 2 outgoing")
    in1, in2 = incoming
    out1, out2 = outgoing
    legs = (in1, in2, out1, out2)
    channels = LoopChannel[]

    # Identify the two fermion lines by matching same-flavor pairs.
    # For ee→μμ: line 1 = (in1=e⁻, in2=e⁺), line 2 = (out1=μ⁻, out2=μ⁺)
    # The incoming pair shares one flavor, the outgoing pair shares another.
    line1 = (in1, in2)   # electron line (particle + antiparticle)
    line2 = (out1, out2) # muon line (particle + antiparticle)

    # Each vertex connects one fermion from line1 + one boson
    # (and similarly for line2). Check that the vertices exist.
    for boson in boson_fields(model)
        bname = boson.name
        # Line 1 vertex: both legs must form valid vertices with the boson
        v_line1 = _has_vertex(rules, in1.field_name, in1.field_name, bname)
        # Line 2 vertex: both legs must form valid vertices with the boson
        v_line2 = _has_vertex(rules, out1.field_name, out1.field_name, bname)
        v_line1 && v_line2 || continue

        # Internal fields: boson₁, line1_fermion, boson₂, line2_fermion
        internals = (bname, in1.field_name, bname, out1.field_name)

        # Direct box: photon connects in1(e⁻) to out1(μ⁻) side
        push!(channels, LoopChannel(:direct_box, internals, legs, :q))
        # Crossed box: photon connects in1(e⁻) to out2(μ⁺) side
        push!(channels, LoopChannel(:crossed_box, internals, legs, :q))
    end
    channels
end
