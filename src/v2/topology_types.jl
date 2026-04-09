# Feynman diagram topology type definition.
#
# A FeynmanTopology is an abstract graph: adjacency matrix + vertex degrees.
# External vertices have degree 1 (indices 1..n_ext).
# Internal vertices have degree ≥ 3 (indices n_ext+1..n).

using LinearAlgebra: diag

struct FeynmanTopology
    n_ext::Int
    adj::Matrix{Int8}
    vdeg::Vector{Int}
end

n_vertices(t::FeynmanTopology) = size(t.adj, 1)
n_internal(t::FeynmanTopology) = n_vertices(t) - t.n_ext

function n_propagators(t::FeynmanTopology)
    n = n_vertices(t)
    total = 0
    for i in 1:n, j in i:n
        total += t.adj[i, j]
    end
    total
end

n_loops(t::FeynmanTopology) = n_propagators(t) - n_vertices(t) + 1
