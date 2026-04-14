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
    count_diagrams(model, in_fields, out_fields; loops=0, onepi=false,
                   nosbridge=false, notadpole=false, onshell=false,
                   nosnail=false, onevi=false,
                   noselfloop=false, nodiloop=false, noparallel=false)

Count connected Feynman diagrams for the given process at the specified
loop order. Phase 17c (Session 24): delegates to `QgrafPort.count_diagrams_qg21`,
the Strategy C qgraf port. The legacy implementation (fill-all-entries
+ is_canonical_feynman) is preserved via `_count_diagrams_legacy` for
regression testing.

Filter kwargs mirror qgraf's dflag options. Currently unsupported:
`nosigma`, `cycli`, `onshellx`, `floop`, `bipart`.
"""
function count_diagrams(model::AbstractModel, in_fields::Vector{Symbol},
                        out_fields::Vector{Symbol};
                        loops::Int=0,
                        onepi::Bool=false,
                        nosbridge::Bool=false,
                        notadpole::Bool=false,
                        onshell::Bool=false,
                        nosnail::Bool=false,
                        onevi::Bool=false,
                        noselfloop::Bool=false,
                        nodiloop::Bool=false,
                        noparallel::Bool=false)
    QgrafPort.count_diagrams_qg21(model, in_fields, out_fields;
                                    loops, onepi, nosbridge, notadpole,
                                    onshell, nosnail, onevi, noselfloop,
                                    nodiloop, noparallel)
end

"""
    _count_diagrams_legacy(model, in_fields, out_fields; loops=0, onepi=false)

Pre-Phase-17c implementation of `count_diagrams`: fill-all-entries topology
enumeration + per-topology field assignment counting. Preserved for
regression testing against the new qg21 path.
"""
function _count_diagrams_legacy(model::AbstractModel, in_fields::Vector{Symbol},
                                out_fields::Vector{Symbol}; loops::Int=0, onepi::Bool=false)
    rules = feynman_rules(model)
    n_ext = length(in_fields) + length(out_fields)
    ext_fields_raw = vcat(in_fields, out_fields)

    exp = _expand_model_for_diagen(model)
    ext_fields = _expand_external_fields(ext_fields_raw, exp)

    vertex_degrees = _model_vertex_degrees(rules)
    isempty(vertex_degrees) && return 0

    count = 0
    for partition in _degree_partitions(n_ext, loops, vertex_degrees)
        for topo in _enumerate_topologies(n_ext, partition; onepi)
            count += _count_field_assignments_expanded(topo, ext_fields, exp)
        end
    end
    count
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

# ---- Expanded model for fermion field splitting ----
#
# qgraf uses separate field names for particle/antiparticle (link array).
# We mirror this: :e → :e (particle) + :e_bar (antiparticle).
# Eliminates orientation ambiguity and overcounting correction.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 1.2 (link array)
# "link(i) = i for self-conjugate; link(particle) = antiparticle"

struct ExpandedModel
    all_fields::Vector{Symbol}
    conjugate::Dict{Symbol, Symbol}
    self_conjugate::Set{Symbol}
    vertex_rules::Set{Vector{Symbol}}
end

function _expand_model_for_diagen(model::AbstractModel)
    all_fields = Symbol[]
    conjugate = Dict{Symbol, Symbol}()
    self_conj = Set{Symbol}()

    for f in model_fields(model)
        if f.self_conjugate
            push!(all_fields, f.name)
            conjugate[f.name] = f.name
            push!(self_conj, f.name)
        else
            anti = Symbol(f.name, :_bar)
            push!(all_fields, f.name)
            push!(all_fields, anti)
            conjugate[f.name] = anti
            conjugate[anti] = f.name
        end
    end

    rules = feynman_rules(model)
    vrules = Set{Vector{Symbol}}()
    for (fields, _) in rules.vertices
        expanded = _expand_vertex(collect(fields), model, conjugate)
        push!(vrules, sort(expanded))
    end

    ExpandedModel(all_fields, conjugate, self_conj, vrules)
end

# Expand one vertex rule: each fermion pair (f,f) → (f, f_bar).
# Ref: refs/qgraf/ALGORITHM.md, Section 1.3 (vertex representation)
function _expand_vertex(fields::Vector{Symbol}, model::AbstractModel,
                        conjugate::Dict{Symbol, Symbol})
    result = Symbol[]
    counts = Dict{Symbol, Int}()
    for f in fields
        counts[f] = get(counts, f, 0) + 1
    end
    for (fname, n) in counts
        fobj = get_field(model, fname)
        if fobj.self_conjugate
            append!(result, fill(fname, n))
        else
            # Fermion number conservation: n/2 particle + n/2 anti
            np = n ÷ 2
            append!(result, fill(fname, np))
            append!(result, fill(conjugate[fname], n - np))
        end
    end
    result
end

function _expand_external_fields(ext_fields::Vector{Symbol}, exp::ExpandedModel)
    # Validate: non-self-conjugate fields must appear in even counts
    for f in ext_fields
        f in exp.self_conjugate && continue
        n = count(==(f), ext_fields)
        iseven(n) || error("Non-self-conjugate field :$f appears $n times (must be even)")
    end
    result = similar(ext_fields)
    seen = Dict{Symbol, Int}()
    for (i, f) in enumerate(ext_fields)
        if f in exp.self_conjugate
            result[i] = f
        else
            c = get(seen, f, 0)
            result[i] = isodd(c) ? exp.conjugate[f] : f
            seen[f] = c + 1
        end
    end
    result
end
