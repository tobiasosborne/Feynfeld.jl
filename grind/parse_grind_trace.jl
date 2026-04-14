#  Parse grind_phi3_2L.txt (or any qgraf-grind trace):
#  extract per-emission topology signature (vmap rows) + ps1, group by
#  topology, print per-topology emission count and a topology summary.
#
#  Usage:  julia parse_grind_trace.jl grind_phi3_2L.txt

length(ARGS) == 1 || error("usage: julia parse_grind_trace.jl <trace_file>")
trace_path = ARGS[1]

# Emission record extracted from the trace.
mutable struct Emit
    id::Int
    n::Int
    n_ext::Int
    ps1::Vector{Int}
    vdeg::Vector{Int}
    vmap::Vector{Vector{Int}}   # vmap[v] = neighbour list
end

emissions = Emit[]

let cur = nothing, ids = 0
  for (lineno, line) in enumerate(eachline(trace_path))
    if startswith(line, "GRIND ===== EMIT_DIAGRAM")
        if cur !== nothing; push!(emissions, cur); end
        ids += 1
        cur = Emit(ids, 0, 0, Int[], Int[], Vector{Vector{Int}}())
    elseif startswith(line, "GRIND EMIT n=")
        m = match(r"n=(\d+) n_ext=(\d+)", line)
        cur.n     = parse(Int, m.captures[1])
        cur.n_ext = parse(Int, m.captures[2])
        resize!(cur.vmap, cur.n)
        for i in 1:cur.n; cur.vmap[i] = Int[]; end
    elseif startswith(line, "GRIND EMIT ps1=")
        s = split(split(line, '=', limit=2)[2])
        cur.ps1 = parse.(Int, s)
    elseif startswith(line, "GRIND EMIT vdeg=")
        s = split(split(line, '=', limit=2)[2])
        cur.vdeg = parse.(Int, s)
    elseif (m = match(r"^GRIND EMIT vmap\[(\d+),1\.\.vdeg\]=(.*)$", line)) !== nothing
        v = parse(Int, m.captures[1])
        nbs = parse.(Int, split(m.captures[2]))
        cur.vmap[v] = nbs
    end
  end
  cur !== nothing && push!(emissions, cur)
end

println("parsed $(length(emissions)) emissions")

# Build topology signature: a tuple-of-tuples of vmap, ignoring ps1 (which is
# ext-labelling, not graph structure).  For qgraf this groups emissions by
# underlying canonical topology.
function topo_sig(e::Emit)
    # qgraf already canonicalises topology so same-topology emissions have
    # identical vmap rows.  Use the full vmap as the signature.
    (e.vdeg, tuple(map(tuple, e.vmap)...))
end

# Also extract an adjacency matrix (gam) from vmap for readability.
function gam_from_vmap(e::Emit)
    g = zeros(Int, e.n, e.n)
    for v in 1:e.n, nb in e.vmap[v]
        g[v, nb] += 1
    end
    # symmetric by construction; divide by 2 on off-diagonal — but we store
    # raw count; gam[i,j] off-diag is count, gam[i,i] is self-loop half-edges.
    # Convert to pseudograph convention: g[i,i] /= 2.
    for i in 1:e.n
        g[i, i] = g[i, i] ÷ 2
    end
    g
end

# Classify edge types per topology for bug diagnosis.
function topo_features(e::Emit)
    g = gam_from_vmap(e)
    n = e.n; ne = e.n_ext
    n_self = 0         # self-loops on internal vertices
    n_parallel = 0     # pairs (i,j) i<j with gam > 1 (parallel edges)
    n_di = 0           # pairs (i,j) i<j with gam == 2 (exactly double edge)
    n_tri = 0          # pairs (i,j) i<j with gam == 3 (triple)
    for i in (ne+1):n
        n_self += g[i, i]
    end
    for i in 1:n, j in (i+1):n
        if g[i, j] >= 2
            n_parallel += 1
            n_di  += (g[i, j] == 2)
            n_tri += (g[i, j] == 3)
        end
    end
    (; n_self, n_parallel, n_di, n_tri)
end

# Bucket by topology signature
buckets = Dict{Any, Vector{Emit}}()
for e in emissions
    push!(get!(Vector{Emit}, buckets, topo_sig(e)), e)
end

println("\n=== $(length(buckets)) distinct topologies ===")
topo_list = collect(pairs(buckets))
sort!(topo_list; by = p -> -length(p.second))

for (i, (sig, emits)) in enumerate(topo_list)
    e = emits[1]
    g = gam_from_vmap(e)
    f = topo_features(e)
    println("\nTOPO #$i  #emissions=$(length(emits))  features=$(f)")
    println("  vdeg = $(e.vdeg)")
    # print upper-triangle gam with non-zero entries
    edges = String[]
    for i1 in 1:e.n, j1 in i1:e.n
        if g[i1, j1] > 0
            push!(edges, "($(i1),$(j1))=$(g[i1,j1])")
        end
    end
    println("  gam = $(join(edges, " "))")
    # all distinct ps1 values in this bucket
    ps1_set = unique(map(x -> x.ps1, emits))
    println("  distinct ps1 count = $(length(ps1_set)); e.g. $(ps1_set[1])")
end

println("\nTotal: $(sum(length(v) for v in values(buckets))) emissions across $(length(buckets)) topologies")
