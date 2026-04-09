# Vertex-degree partition enumeration for Feynman diagram generation.
#
# Given n_ext external legs, L loops, and allowed vertex degrees,
# enumerate all valid vertex-degree sequences ρ(k) such that:
#   Σ_k k·ρ(k) + n_ext = 2·P  (total half-edges = 2×propagators + external)
#   P - V + 1 = L              (Euler formula for loops)
# where V = Σ_k ρ(k) and P = (Σ_k k·ρ(k) + n_ext) / 2.
#
# Ref: refs/qgraf/ALGORITHM.md, Section 2

"""
    DegreePartition

A vertex-degree partition: how many internal vertices of each degree.
"""
struct DegreePartition
    # degree => count, e.g. Dict(3 => 4) means 4 vertices of degree 3
    counts::Dict{Int, Int}
end

n_vertices(dp::DegreePartition) = sum(values(dp.counts); init=0)
total_half_edges(dp::DegreePartition) = sum(d * c for (d, c) in dp.counts; init=0)

"""
    _model_vertex_degrees(rules) → Set{Int}

Extract the set of vertex degrees available in the model.
"""
function _model_vertex_degrees(rules::FeynmanRules)
    Set(length(k) for k in keys(rules.vertices))
end

"""
    _degree_partitions(n_ext, loops, vertex_degrees) → Vector{DegreePartition}

Enumerate all valid vertex-degree partitions for n_ext external legs
at the given loop order, using only the specified vertex degrees.

Constraint: Σ_k (k-2)·ρ(k) = n_ext + 2·(loops - 1)
"""
function _degree_partitions(n_ext::Int, loops::Int, vertex_degrees::Set{Int})
    target = n_ext + 2 * (loops - 1)
    # target = Σ_k (k-2)·ρ(k)
    # For degree 3: each vertex contributes 1
    # For degree 4: each vertex contributes 2
    # etc.

    degs = sort(collect(vertex_degrees))
    result = DegreePartition[]
    _partition_recurse!(result, degs, 1, target, Dict{Int,Int}(), n_ext, loops)
    result
end

function _partition_recurse!(result, degs, idx, remaining, counts, n_ext, loops)
    if remaining == 0
        # Check that the partition gives the correct loop count
        dp = DegreePartition(copy(counts))
        n_v_internal = n_vertices(dp)
        n_v_total = n_v_internal + n_ext  # ALL vertices: internal + external
        total_he = total_half_edges(dp) + n_ext
        total_he % 2 == 0 || return  # must be even
        n_prop = total_he ÷ 2
        actual_loops = n_prop - n_v_total + 1  # Euler: L = P - V + 1
        actual_loops == loops || return
        push!(result, dp)
        return
    end
    remaining < 0 && return
    idx > length(degs) && return

    d = degs[idx]
    contrib = d - 2  # each vertex of degree d contributes d-2 to the target
    contrib <= 0 && return  # degree ≤ 2 can't contribute positively

    # Try 0, 1, 2, ... vertices of degree d
    max_count = remaining ÷ contrib
    for c in 0:max_count
        if c > 0
            counts[d] = c
        end
        _partition_recurse!(result, degs, idx + 1, remaining - c * contrib, counts, n_ext, loops)
        if c > 0
            delete!(counts, d)
        end
    end
end
