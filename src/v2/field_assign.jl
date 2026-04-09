# Field assignment: map abstract topology to concrete Feynman diagrams.
#
# Given a topology and a model, count valid field assignments.
# External edges are fixed. Internal edges are assigned by backtracking.
# Fermion flow (particle/antiparticle) is tracked at each vertex.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 4

"""
    _count_field_assignments(topo, ext_fields, is_anti, rules, model) → Int

Count valid field assignments. External edges fixed, internal variable.
Fermion orientations tracked. Closed fermion loop overcounting corrected.
"""
function _count_field_assignments(topo::FeynmanTopology, ext_fields::Vector{Symbol},
                                  is_anti::Vector{Bool},
                                  rules::FeynmanRules, model::AbstractModel)
    n = n_vertices(topo)
    n_ext = topo.n_ext

    edges = Tuple{Int,Int}[]
    edge_fixed = Symbol[]
    edge_anti_at_i = Union{Bool, Nothing}[]
    for i in 1:n, j in i:n
        for _ in 1:topo.adj[i, j]
            push!(edges, (i, j))
            if i <= n_ext
                push!(edge_fixed, ext_fields[i])
                push!(edge_anti_at_i, is_anti[i])
            else
                push!(edge_fixed, :none)
                push!(edge_anti_at_i, nothing)
            end
        end
    end

    all_fields = [f.name for f in model_fields(model)]
    vertex_rules = _vertex_rule_multisets(rules)

    edge_fields = copy(edge_fixed)
    count = Ref(0//1)  # Rational for fermion loop correction
    _assign_edges!(count, edge_fields, edge_fixed, edge_anti_at_i, edges, 1,
                   topo, ext_fields, is_anti, all_fields, vertex_rules, model)
    Int(count[])
end

"""
    _vertex_rule_multisets(rules) → Set{Vector{Symbol}}
"""
function _vertex_rule_multisets(rules::FeynmanRules)
    result = Set{Vector{Symbol}}()
    for (fields, _) in rules.vertices
        push!(result, sort(collect(fields)))
    end
    result
end

function _assign_edges!(count, edge_fields, edge_fixed, edge_anti_at_i, edges,
                        idx, topo, ext_fields, is_anti, all_fields,
                        vertex_rules, model)
    n_edges = length(edges)
    n_ext = topo.n_ext
    n = n_vertices(topo)

    if idx > n_edges
        for v in (n_ext + 1):n
            _check_vertex_with_flow(v, topo, edge_fields, edge_anti_at_i, edges,
                                    ext_fields, is_anti, vertex_rules, model) || return
        end
        # Correct for closed fermion loop overcounting
        n_fl = _count_closed_fermion_loops(topo, edge_fields, edge_anti_at_i,
                                            edges, is_anti, model)
        count[] += 1 // (1 << n_fl)
        return
    end

    if edge_fixed[idx] != :none
        _assign_edges!(count, edge_fields, edge_fixed, edge_anti_at_i, edges,
                        idx + 1, topo, ext_fields, is_anti, all_fields,
                        vertex_rules, model)
        return
    end

    i, j = edges[idx]
    for fname in all_fields
        field_obj = get_field(model, fname)
        if field_obj.self_conjugate
            edge_fields[idx] = fname
            edge_anti_at_i[idx] = false
            _try_assignment!(count, edge_fields, edge_fixed, edge_anti_at_i,
                             edges, idx, i, j, n_ext, topo, ext_fields, is_anti,
                             all_fields, vertex_rules, model)
        else
            for anti_at_i in (false, true)
                edge_fields[idx] = fname
                edge_anti_at_i[idx] = anti_at_i
                _try_assignment!(count, edge_fields, edge_fixed, edge_anti_at_i,
                                 edges, idx, i, j, n_ext, topo, ext_fields,
                                 is_anti, all_fields, vertex_rules, model)
            end
        end
    end
    edge_fields[idx] = :none
    edge_anti_at_i[idx] = nothing
end

function _try_assignment!(count, edge_fields, edge_fixed, edge_anti_at_i,
                          edges, idx, i, j, n_ext, topo, ext_fields, is_anti,
                          all_fields, vertex_rules, model)
    valid = true
    for v in (i, j)
        v > n_ext || continue
        _partial_vertex_ok_flow(v, idx, topo, edge_fields, edge_anti_at_i,
                                edges, ext_fields, is_anti, vertex_rules,
                                model) || (valid = false; break)
    end
    valid && _assign_edges!(count, edge_fields, edge_fixed, edge_anti_at_i,
                            edges, idx + 1, topo, ext_fields, is_anti,
                            all_fields, vertex_rules, model)
end

# Vertex checking and fermion flow functions are in vertex_check.jl
# Fermion loop counting is in fermion_flow.jl
