# Fermion flow: count closed fermion loops for overcounting correction.
#
# Each closed fermion loop can be traversed in 2 directions, giving
# 2× overcounting per loop. We divide by 2^n_closed_loops.
#
# A closed loop = connected component of internal fermion edges
# that does NOT connect to any external fermion leg.

"""
    _count_closed_fermion_loops(topo, edge_fields, edge_anti_at_i, edges,
                                is_anti, model) → Int

Count closed fermion loops (internal fermion cycles with no external legs).
"""
function _count_closed_fermion_loops(topo, edge_fields, edge_anti_at_i, edges,
                                     is_anti, model)
    n_ext = topo.n_ext

    # Internal fermion edges
    fermion_edges = Int[]
    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        i > n_ext && j > n_ext || continue
        field_obj = get_field(model, edge_fields[idx])
        field_obj.self_conjugate && continue
        push!(fermion_edges, idx)
    end
    isempty(fermion_edges) && return 0

    # Union-find on shared internal vertices
    parent = Dict(idx => idx for idx in fermion_edges)
    function uf_find(x)
        while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x
    end

    v_to_edges = Dict{Int, Vector{Int}}()
    for idx in fermion_edges
        i, j = edges[idx]
        push!(get!(Vector{Int}, v_to_edges, i), idx)
        if i != j
            push!(get!(Vector{Int}, v_to_edges, j), idx)
        end
    end
    for (_, vedges) in v_to_edges
        for k in 2:length(vedges)
            ra, rb = uf_find(vedges[1]), uf_find(vedges[k])
            ra != rb && (parent[ra] = rb)
        end
    end

    # Which internal vertices are reached by external fermion legs?
    ext_fermion_int_verts = Set{Int}()
    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        if i <= n_ext
            field_obj = get_field(model, edge_fields[idx])
            !field_obj.self_conjugate && push!(ext_fermion_int_verts, j)
        end
    end

    # Group and check each component
    comp_edges = Dict{Int, Vector{Int}}()
    for idx in fermion_edges
        root = uf_find(idx)
        push!(get!(Vector{Int}, comp_edges, root), idx)
    end

    n_closed = 0
    for (_, edge_list) in comp_edges
        # Does this component touch external fermion vertices?
        touches_ext = false
        for idx in edge_list
            i, j = edges[idx]
            if i in ext_fermion_int_verts || j in ext_fermion_int_verts
                touches_ext = true
                break
            end
        end
        touches_ext || (n_closed += 1)
    end
    n_closed
end
