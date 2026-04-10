# Field assignment: map abstract topology to concrete Feynman diagrams.
#
# Given a topology and an expanded model, count valid field assignments.
# External edges are fixed. Internal edges are assigned by backtracking.
#
# With expanded fields (particle ≠ antiparticle), each field name carries
# full identity. Vertex i sees field f, vertex j sees conjugate(f).
# Canonical ordering on parallel edges prevents multi-edge permutation.
# Closed fermion loop overcounting (self-loops) corrected via ÷2^n.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 4

function _count_field_assignments_expanded(topo::FeynmanTopology,
        ext_fields::Vector{Symbol}, exp)
    n = n_vertices(topo)
    n_ext = topo.n_ext

    edges = Tuple{Int,Int}[]
    edge_fixed = Symbol[]
    for i in 1:n, j in i:n
        for _ in 1:topo.adj[i, j]
            push!(edges, (i, j))
            push!(edge_fixed, i <= n_ext ? ext_fields[i] : :none)
        end
    end

    edge_fields = copy(edge_fixed)
    count = Ref(0//1)
    _assign_expanded!(count, edge_fields, edge_fixed, edges, 1,
                      topo, exp.all_fields, exp.vertex_rules,
                      exp.conjugate, exp.self_conjugate)
    Int(count[])
end

function _assign_expanded!(count, edge_fields, edge_fixed, edges, idx,
                           topo, all_fields, vertex_rules, conjugate, self_conj)
    if idx > length(edges)
        n_ext = topo.n_ext
        for v in (n_ext + 1):n_vertices(topo)
            _check_vertex_expanded(v, topo, edge_fields, edges,
                                   vertex_rules, conjugate) || return
        end
        # Closed fermion loops are counted in both orientations → divide by 2 each
        n_fl = _count_closed_loops_expanded(topo, edge_fields, edges, self_conj)
        count[] += 1 // (1 << n_fl)
        return
    end

    if edge_fixed[idx] != :none
        _assign_expanded!(count, edge_fields, edge_fixed, edges,
                          idx + 1, topo, all_fields, vertex_rules,
                          conjugate, self_conj)
        return
    end

    i, j = edges[idx]
    n_ext = topo.n_ext
    # Canonical ordering on parallel edges: if this edge has the same
    # (i,j) as the previous one, only try fields ≥ the previous field.
    # Prevents counting permutations of indistinguishable multi-edges.
    min_field = :_  # underscore sorts before all letters
    if idx > 1 && edges[idx] == edges[idx - 1] && edge_fixed[idx] == :none &&
                   edge_fixed[idx - 1] == :none
        min_field = edge_fields[idx - 1]
    end
    for fname in all_fields
        fname >= min_field || continue
        edge_fields[idx] = fname
        valid = true
        for v in (i, j)
            v > n_ext || continue
            _partial_vertex_ok_expanded(v, idx, topo, edge_fields, edges,
                                        vertex_rules, conjugate) || (valid = false; break)
        end
        if valid
            _assign_expanded!(count, edge_fields, edge_fixed, edges,
                              idx + 1, topo, all_fields, vertex_rules,
                              conjugate, self_conj)
        end
    end
    edge_fields[idx] = :none
end

# Count closed fermion loops: connected components of non-self-conjugate
# internal edges whose vertices don't touch any external vertex.
# Simpler than the old fermion_flow.jl because expanded fields make
# fermion detection trivial (just check self_conj membership).
function _count_closed_loops_expanded(topo, edge_fields, edges, self_conj)
    n_ext = topo.n_ext

    # Internal non-self-conjugate edges
    fedges = Int[]
    for (idx, (i, j)) in enumerate(edges)
        edge_fields[idx] == :none && continue
        i > n_ext && j > n_ext || continue
        edge_fields[idx] in self_conj && continue
        push!(fedges, idx)
    end
    isempty(fedges) && return 0

    # Which internal vertices carry an external FERMION leg?
    # (Photon legs don't break fermion loop closure)
    ext_adj = Set{Int}()
    for (idx, (i, j)) in enumerate(edges)
        i <= n_ext || continue
        edge_fields[idx] in self_conj && continue
        push!(ext_adj, j)
    end

    # Union-find on shared internal vertices
    par = Dict(idx => idx for idx in fedges)
    uf(x) = (while par[x] != x; par[x] = par[par[x]]; x = par[x]; end; x)

    vmap = Dict{Int, Vector{Int}}()
    for idx in fedges
        i, j = edges[idx]
        push!(get!(Vector{Int}, vmap, i), idx)
        i != j && push!(get!(Vector{Int}, vmap, j), idx)
    end
    for (_, ve) in vmap
        for k in 2:length(ve)
            ra, rb = uf(ve[1]), uf(ve[k])
            ra != rb && (par[ra] = rb)
        end
    end

    # Count components not touching external-adjacent vertices.
    # Skip multi-edge components (all edges between same distinct pair) —
    # those are already handled by the canonical edge ordering.
    comps = Dict{Int, Vector{Int}}()
    for idx in fedges
        push!(get!(Vector{Int}, comps, uf(idx)), idx)
    end
    n_closed = 0
    for (_, elist) in comps
        touches = false
        for idx in elist
            i, j = edges[idx]
            (i in ext_adj || j in ext_adj) && (touches = true; break)
        end
        touches && continue
        # Check if this is a multi-edge (all edges between same distinct pair)
        if length(elist) >= 2
            i0, j0 = edges[elist[1]]
            if i0 != j0 && all(edges[idx] == (i0, j0) for idx in elist)
                continue  # multi-edge: canonical ordering handles it
            end
        end
        n_closed += 1
    end
    n_closed
end
