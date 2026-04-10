# Vertex constraint checking for expanded-field assignment.
#
# With expanded field names (particle ≠ antiparticle), vertex checking
# is a sorted-multiset comparison against expanded vertex rules.
# Key: vertex i sees field f, vertex j sees conjugate(f).
#
# Ref: refs/qgraf/ALGORITHM.md, Section 4.1
# "pmap(vmap(i1,1), lmap(i1,1)) = link(0)+j1" — j-end sees conjugate

function _collect_fields_expanded(v, topo, edge_fields, edges,
                                  conjugate::Dict{Symbol, Symbol})
    result = Symbol[]
    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        f = edge_fields[idx]
        if i == j && i == v
            push!(result, f)
            push!(result, conjugate[f])
        elseif i == v
            push!(result, f)
        elseif j == v
            push!(result, conjugate[f])
        end
    end
    result
end

function _check_vertex_expanded(v, topo, edge_fields, edges,
                                vertex_rules, conjugate)
    fields_at_v = _collect_fields_expanded(v, topo, edge_fields, edges, conjugate)
    sort(fields_at_v) in vertex_rules
end

function _partial_vertex_ok_expanded(v, up_to_edge, topo, edge_fields, edges,
                                     vertex_rules, conjugate)
    assigned = Symbol[]
    for (idx, (i, j)) in enumerate(edges)
        (i == v || j == v) || continue
        idx <= up_to_edge || continue
        edge_fields[idx] == :none && continue
        f = edge_fields[idx]
        if i == j && i == v
            push!(assigned, f)
            push!(assigned, conjugate[f])
        elseif i == v
            push!(assigned, f)
        elseif j == v
            push!(assigned, conjugate[f])
        end
    end
    isempty(assigned) && return true
    length(assigned) == topo.vdeg[v] && return sort(assigned) in vertex_rules
    sorted = sort(assigned)
    for rule in vertex_rules
        length(rule) == topo.vdeg[v] || continue
        _is_submultiset(sorted, rule) && return true
    end
    false
end

# Check if `sub` (sorted) is a sub-multiset of `sup` (sorted).
function _is_submultiset(sub::Vector{Symbol}, sup::Vector{Symbol})
    j = 1
    for s in sub
        while j <= length(sup) && sup[j] < s
            j += 1
        end
        j > length(sup) && return false
        sup[j] == s || return false
        j += 1
    end
    true
end
