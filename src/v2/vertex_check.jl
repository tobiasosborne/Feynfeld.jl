# Vertex constraint checking for field assignment.
#
# At each vertex, the multiset of field names on its half-edges must
# match a vertex rule in the model.

"""
    _check_vertex_fields(v, topo, edge_fields, edges, vertex_rules) → Bool

Check if vertex v's field multiset matches a vertex rule (no orientation check).
"""
function _check_vertex_fields(v, topo, edge_fields, edges, vertex_rules)
    fields_at_v = _collect_fields_at_vertex(v, topo, edge_fields, edges)
    sort(fields_at_v) in vertex_rules
end

"""
    _partial_vertex_fields_ok(v, up_to_edge, topo, edge_fields, edges, vertex_rules) → Bool

Partial check: assigned fields so far compatible with some vertex rule?
"""
function _partial_vertex_fields_ok(v, up_to_edge, topo, edge_fields, edges, vertex_rules)
    assigned = Symbol[]
    for (idx, (i, j)) in enumerate(edges)
        (i == v || j == v) || continue
        idx <= up_to_edge || continue
        edge_fields[idx] == :none && continue
        if i == j && i == v
            push!(assigned, edge_fields[idx])
            push!(assigned, edge_fields[idx])
        else
            push!(assigned, edge_fields[idx])
        end
    end
    isempty(assigned) && return true
    length(assigned) == topo.vdeg[v] && return sort(assigned) in vertex_rules
    sorted_assigned = sort(assigned)
    for rule in vertex_rules
        length(rule) == topo.vdeg[v] || continue
        _is_submultiset(sorted_assigned, rule) && return true
    end
    false
end

function _collect_fields_at_vertex(v, topo, edge_fields, edges)
    result = Symbol[]
    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        if i == j && i == v
            push!(result, edge_fields[idx])
            push!(result, edge_fields[idx])
        elseif i == v || j == v
            push!(result, edge_fields[idx])
        end
    end
    result
end

"""
    _check_vertex_with_flow(v, ...) → Bool

Check vertex v: field multiset matches a vertex rule AND fermion flow valid.
"""
function _check_vertex_with_flow(v, topo, edge_fields, edge_anti_at_i, edges,
                                  ext_fields, is_anti, vertex_rules, model)
    fields_at_v = Symbol[]
    n_particle = 0
    n_anti = 0
    n_ext = topo.n_ext

    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        if i == j && i == v
            f = edge_fields[idx]
            push!(fields_at_v, f)
            push!(fields_at_v, f)
            field_obj = get_field(model, f)
            if !field_obj.self_conjugate
                n_particle += 1
                n_anti += 1
            end
        elseif i == v || j == v
            f = edge_fields[idx]
            push!(fields_at_v, f)
            field_obj = get_field(model, f)
            if !field_obj.self_conjugate
                anti_here = _is_anti_at_vertex(v, i, j, idx, edge_anti_at_i,
                                                ext_fields, is_anti, n_ext)
                anti_here ? (n_anti += 1) : (n_particle += 1)
            end
        end
    end

    sort(fields_at_v) in vertex_rules || return false
    n_particle == n_anti || return false
    true
end

"""
    _partial_vertex_ok_flow(v, up_to_edge, ...) → Bool

Partial check: are the assigned edges so far compatible with some vertex rule
and consistent fermion flow?
"""
function _partial_vertex_ok_flow(v, up_to_edge, topo, edge_fields, edge_anti_at_i,
                                  edges, ext_fields, is_anti, vertex_rules, model)
    assigned = Symbol[]
    n_particle = 0
    n_anti = 0
    n_ext = topo.n_ext

    for (idx, (i, j)) in enumerate(edges)
        (i == v || j == v) || continue
        idx <= up_to_edge || continue
        edge_fields[idx] == :none && continue
        f = edge_fields[idx]
        if i == j && i == v
            push!(assigned, f)
            push!(assigned, f)
            field_obj = get_field(model, f)
            if !field_obj.self_conjugate
                n_particle += 1
                n_anti += 1
            end
        else
            push!(assigned, f)
            field_obj = get_field(model, f)
            if !field_obj.self_conjugate
                anti_here = _is_anti_at_vertex(v, i, j, idx, edge_anti_at_i,
                                                ext_fields, is_anti, n_ext)
                anti_here ? (n_anti += 1) : (n_particle += 1)
            end
        end
    end

    isempty(assigned) && return true
    n_remaining = topo.vdeg[v] - length(assigned)
    abs(n_particle - n_anti) > n_remaining && return false
    length(assigned) == topo.vdeg[v] && return sort(assigned) in vertex_rules

    sorted_assigned = sort(assigned)
    for rule in vertex_rules
        length(rule) == topo.vdeg[v] || continue
        _is_submultiset(sorted_assigned, rule) && return true
    end
    false
end

"""
    _is_anti_at_vertex(v, i, j, idx, ...) → Bool

Is vertex v seeing the antiparticle end of edge idx?
"""
function _is_anti_at_vertex(v, i, j, idx, edge_anti_at_i, ext_fields, is_anti, n_ext)
    # External edge: both ends see the same particle/antiparticle identity
    i <= n_ext && return is_anti[i]
    j <= n_ext && return is_anti[j]
    # Internal edge: vertex i sees edge_anti_at_i, vertex j sees opposite
    aai = edge_anti_at_i[idx]
    aai === nothing && return false
    (v == i) ? aai : !aai
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
