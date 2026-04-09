# Layer 3: Algorithmic Feynman diagram generation.
#
# Replaces the hard-coded channel enumeration (channels.jl, loop_channels.jl)
# with an algorithmic approach based on qgraf's adjacency-matrix method.
#
# Algorithm: for each valid vertex-degree partition, enumerate all
# non-isomorphic topologies via constrained adjacency matrix fill,
# then assign fields from the model. See refs/qgraf/ALGORITHM.md.
#
# Ref: P. Nogueira, J. Comput. Phys. 105 (1993) 279
# "Automatic Feynman graph generation"

using Combinatorics: partitions, multiset_permutations

"""
    count_diagrams(model, in_fields, out_fields; loops=0, onepi=false)

Count connected Feynman diagrams for the given process at the specified
loop order. Field names are Symbols matching the model's field names.

For fermions, the same name is used for particle and antiparticle.
The generator infers particle/antiparticle from the standard convention:
  - incoming: first leg = particle, second = antiparticle
  - outgoing: first leg = particle, second = antiparticle
External legs are paired as: in1(particle), in2(anti), out1(particle), out2(anti)
for 2→2. For general processes, the `is_anti` vector marks antiparticles.
"""
function count_diagrams(model::AbstractModel, in_fields::Vector{Symbol},
                        out_fields::Vector{Symbol}; loops::Int=0, onepi::Bool=false)
    rules = feynman_rules(model)
    n_ext = length(in_fields) + length(out_fields)
    ext_fields = vcat(in_fields, out_fields)

    # Build is_anti vector: for each external leg, is it an antiparticle?
    # Convention: for each field that appears as a pair in in/out,
    # the second occurrence is the antiparticle.
    is_anti = _infer_antiparticles(ext_fields, model)

    # Collect vertex degrees available in the model
    vertex_degrees = _model_vertex_degrees(rules)
    isempty(vertex_degrees) && return 0

    # Enumerate all valid topologies and count field-compatible diagrams
    count = 0
    for partition in _degree_partitions(n_ext, loops, vertex_degrees)
        for topo in _enumerate_topologies(n_ext, partition; onepi)
            count += _count_field_assignments(topo, ext_fields, is_anti, rules, model)
        end
    end
    count
end

"""
    _infer_antiparticles(ext_fields, model) → Vector{Bool}

Infer which external legs are antiparticles based on field properties.
For self-conjugate fields (bosons like photon, scalar like phi): false.
For non-self-conjugate fields (fermions): even-numbered occurrences
of the same field are antiparticles (convention: particle, anti, particle, anti).
"""
function _infer_antiparticles(ext_fields::Vector{Symbol}, model::AbstractModel)
    is_anti = fill(false, length(ext_fields))
    seen_count = Dict{Symbol, Int}()
    for (i, f) in enumerate(ext_fields)
        field_obj = get_field(model, f)
        if field_obj.self_conjugate
            is_anti[i] = false
        else
            c = get(seen_count, f, 0)
            is_anti[i] = isodd(c)  # 0th=particle, 1st=anti, 2nd=particle, ...
            seen_count[f] = c + 1
        end
    end
    is_anti
end

"""
    generate_tree_channels(model, rules, incoming, outgoing)

Generate tree-level channels using the algorithmic diagram generator.
Returns Vector{TreeChannel} compatible with existing pipeline.
"""
function generate_tree_channels(model::AbstractModel, rules::FeynmanRules,
                                incoming::Vector{ExternalLeg},
                                outgoing::Vector{ExternalLeg})
    # For 2→2 tree-level, delegate to existing tree_channels for now.
    # TODO: replace with full algorithmic generation once tested.
    tree_channels(model, rules, incoming, outgoing)
end
