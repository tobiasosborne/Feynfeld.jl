#  Compare topology inventories between qgraf (grind_phi3_2L.txt) and
#  Julia (julia_phi3_2L.txt).  Strategy: build an isomorphism-invariant
#  signature for each topology (sorted degree sequence of its line graph
#  + list of edge-multiplicities) and cross-tab.
#
#  This is a diagnostic — not perfect isomorphism, but strong enough to
#  highlight the extra Julia topologies.

# ── helpers ──────────────────────────────────────────────────────────
function canon_sig(n::Int, n_ext::Int, gam::Matrix{Int})
    # Iso-invariant signature: for each vertex, (degree, xn, self-loops,
    # sorted multiset of neighbour "labels" where each neighbour-label
    # includes the edge multiplicity).  Then sort all vertex signatures.
    #
    # Uses 2 rounds of Weisfeiler-Lehman-like refinement.
    deg      = [sum(gam[i, :]) + gam[i, i] for i in 1:n]   # half-edge degree
    xn       = [sum(j <= n_ext ? gam[i, j] : 0 for j in 1:n) for i in 1:n]
    xg_diag  = [gam[i, i] for i in 1:n]
    lab0     = [(deg[i], xn[i], xg_diag[i], i <= n_ext ? :ext : :int) for i in 1:n]
    # round 1: include sorted multiset of neighbour labels
    lab1 = Vector{Any}(undef, n)
    for i in 1:n
        nbs = Tuple[]
        for j in 1:n
            mult = j == i ? gam[i, j] ÷ 2 : gam[i, j]
            mult > 0 || continue
            push!(nbs, (mult, lab0[j]))
        end
        sort!(nbs)
        lab1[i] = (lab0[i], tuple(nbs...))
    end
    # round 2: refine once more
    lab2 = Vector{Any}(undef, n)
    for i in 1:n
        nbs = Tuple[]
        for j in 1:n
            mult = j == i ? gam[i, j] ÷ 2 : gam[i, j]
            mult > 0 || continue
            push!(nbs, (mult, lab1[j]))
        end
        sort!(nbs)
        lab2[i] = (lab1[i], tuple(nbs...))
    end
    tuple(sort(lab2)...)
end

# ── parse qgraf trace (half-edge convention for diag) ────────────────
function parse_qgraf(path)
    emissions = []
    let cur = nothing
        for line in eachline(path)
            if startswith(line, "GRIND ===== EMIT_DIAGRAM")
                if cur !== nothing; push!(emissions, cur); end
                cur = Dict{Symbol, Any}(:vmap => Dict{Int, Vector{Int}}())
            elseif (m = match(r"^GRIND EMIT n=(\d+) n_ext=(\d+)", line)) !== nothing
                cur[:n] = parse(Int, m.captures[1]); cur[:n_ext] = parse(Int, m.captures[2])
            elseif startswith(line, "GRIND EMIT vdeg=")
                cur[:vdeg] = parse.(Int, split(split(line, '=', limit=2)[2]))
            elseif startswith(line, "GRIND EMIT ps1=")
                cur[:ps1] = parse.(Int, split(split(line, '=', limit=2)[2]))
            elseif (m = match(r"^GRIND EMIT vmap\[(\d+),1\.\.vdeg\]=(.*)$", line)) !== nothing
                v = parse(Int, m.captures[1])
                cur[:vmap][v] = parse.(Int, split(m.captures[2]))
            end
        end
        cur !== nothing && push!(emissions, cur)
    end
    emissions
end

function qgraf_gam(e)
    n = e[:n]
    g = zeros(Int, n, n)
    for v in 1:n, nb in e[:vmap][v]
        g[v, nb] += 1
    end
    # qgraf vmap counts BOTH half-edges of self-loops, so g[i,i] already
    # in half-edge convention.  off-diagonal g[i,j] = edge multiplicity.
    # For off-diag: vmap[i] contains j `gam(i,j)` times → g[i,j] = gam.
    # And g[j,i] also = gam.  So g is symmetric in off-diag.
    g
end

# ── parse Julia dump (already half-edge conv for xg diag) ────────────
function parse_julia(path)
    topos = []
    let cur = nothing
        for line in eachline(path)
            if (m = match(r"^JULIA TOPO #(\d+)\s+#emissions=(\d+)\s+\|G\|=(\d+)\s+burnside_contrib=(.*)$", line)) !== nothing
                if cur !== nothing; push!(topos, cur); end
                cur = Dict{Symbol, Any}(:id => parse(Int, m.captures[1]),
                                         :nem => parse(Int, m.captures[2]),
                                         :G => parse(Int, m.captures[3]),
                                         :contrib => m.captures[4])
            elseif (m = match(r"^  features=\(n_self=(\d+), n_parallel=(\d+), n_di=(\d+), n_tri=(\d+)\)", line)) !== nothing
                cur[:n_self]     = parse(Int, m.captures[1])
                cur[:n_parallel] = parse(Int, m.captures[2])
                cur[:n_di]       = parse(Int, m.captures[3])
                cur[:n_tri]      = parse(Int, m.captures[4])
            elseif startswith(line, "  gam=")
                cur[:gam_str] = line[7:end]
            end
        end
        cur !== nothing && push!(topos, cur)
    end
    topos
end

function julia_gam(t)
    # Parse the gam_str of form "(i,j)=m (i,j)=m ..."
    # We don't know n from the string alone; just return a sparse list for
    # now and let canon_sig build the Matrix.
    pairs = Tuple{Int,Int,Int}[]
    for tok in split(t[:gam_str])
        m = match(r"^\((\d+),(\d+)\)=(\d+)$", tok)
        m === nothing && continue
        push!(pairs, (parse(Int, m.captures[1]), parse(Int, m.captures[2]), parse(Int, m.captures[3])))
    end
    n = 0
    for p in pairs
        n = max(n, p[1], p[2])
    end
    g = zeros(Int, n, n)
    for p in pairs
        i, j, mu = p
        if i == j
            g[i, i] = mu   # half-edge count (Julia xg already)
        else
            g[i, j] = mu
            g[j, i] = mu
        end
    end
    n, g
end

# ── main ─────────────────────────────────────────────────────────────
q_emissions = parse_qgraf(joinpath(@__DIR__, "grind_phi3_2L.txt"))
println("qgraf parsed $(length(q_emissions)) emissions")

q_topo_sigs = Dict{Any, Int}()
for e in q_emissions
    g = qgraf_gam(e)
    sig = canon_sig(e[:n], e[:n_ext], g)
    q_topo_sigs[sig] = get(q_topo_sigs, sig, 0) + 1
end
println("qgraf canonical topology count = $(length(q_topo_sigs))")
println("qgraf emission counts per topo: sorted descending")
for (i, (sig, c)) in enumerate(sort(collect(q_topo_sigs); by = p -> -p.second))
    i <= 10 || break
    println("  top $i: $c emissions")
end

j_topos = parse_julia(joinpath(@__DIR__, "julia_phi3_2L.txt"))
println("\nJulia parsed $(length(j_topos)) topologies")

j_topo_sigs = Dict{Any, Int}()
j_sig_to_tid = Dict{Any, Vector{Int}}()
for t in j_topos
    n, g = julia_gam(t)
    # Julia reports a topology with gam in xg half-edge conv; needs padding
    # to n=10 since the string may omit high unused indices — rebuild with n=10.
    g_full = zeros(Int, 10, 10)
    g_full[1:n, 1:n] = g
    sig = canon_sig(10, 4, g_full)
    j_topo_sigs[sig] = get(j_topo_sigs, sig, 0) + 1
    push!(get!(Vector{Int}, j_sig_to_tid, sig), t[:id])
end
println("Julia canonical topology count (by iso-sig) = $(length(j_topo_sigs))")

# Count topologies per iso-class in each
q_per_iso = Dict{Any, Int}()  # sig → number of DISTINCT qgraf canonical topos with this sig
q_sig_by_emissions = Dict{Any, Vector{Any}}()  # sig → list of topology representatives
let seen_topos_q = Set{Any}()
    for e in q_emissions
        g = qgraf_gam(e)
        sig = canon_sig(e[:n], e[:n_ext], g)
        matsig = (e[:vdeg], g)  # exact-matrix sig, per-emission
        if !(matsig in seen_topos_q)
            push!(seen_topos_q, matsig)
            q_per_iso[sig] = get(q_per_iso, sig, 0) + 1
        end
    end
end

j_per_iso = Dict{Any, Int}()
for t in j_topos
    n, g = julia_gam(t)
    g_full = zeros(Int, 10, 10); g_full[1:n, 1:n] = g
    sig = canon_sig(10, 4, g_full)
    j_per_iso[sig] = get(j_per_iso, sig, 0) + 1
end

println("\nqgraf topologies per iso-class: $(sort(collect(values(q_per_iso))))")
println("Julia topologies per iso-class: $(sort(collect(values(j_per_iso))))")

# Iso-classes where Julia has MORE topologies than qgraf
println("\nIso-classes where Julia count > qgraf count:")
for sig in keys(j_per_iso)
    jc = j_per_iso[sig]
    qc = get(q_per_iso, sig, 0)
    if jc > qc
        println("  Julia $jc vs qgraf $qc  (excess $(jc - qc))")
        for t in j_topos
            n, g = julia_gam(t)
            g_full = zeros(Int, 10, 10); g_full[1:n, 1:n] = g
            canon_sig(10, 4, g_full) == sig || continue
            println("    julia id=$(t[:id]) nem=$(t[:nem]) |G|=$(t[:G]) contrib=$(t[:contrib]) nself=$(t[:n_self]) npar=$(t[:n_parallel])")
            println("      gam=$(t[:gam_str])")
        end
    end
end

common = intersect(keys(q_topo_sigs), keys(j_topo_sigs))
only_q = setdiff(keys(q_topo_sigs), keys(j_topo_sigs))
only_j = setdiff(keys(j_topo_sigs), keys(q_topo_sigs))
println("\ncommon signatures: $(length(common))")
println("only-in-qgraf: $(length(only_q))")
println("only-in-julia: $(length(only_j))")

if true
    println("\n=== For each Julia-excess iso-class, show qgraf's canonical form ===")
    # For each iso-class where Julia has > qgraf count, show all qgraf
    # exact-matrix topologies with that iso-class.
    # Build a reverse index from iso-sig → qgraf representatives.
    q_sig_to_reps = Dict{Any, Vector{Any}}()
    let seen_topos_q = Set{Any}()
        for e in q_emissions
            g = qgraf_gam(e)
            sig = canon_sig(e[:n], e[:n_ext], g)
            matsig = (e[:vdeg], g)
            if !(matsig in seen_topos_q)
                push!(seen_topos_q, matsig)
                push!(get!(Vector{Any}, q_sig_to_reps, sig), e)
            end
        end
    end
    for sig in keys(j_per_iso)
        jc = j_per_iso[sig]
        qc = get(q_per_iso, sig, 0)
        jc > qc || continue
        println("\n-- Iso-class with Julia $jc / qgraf $qc --")
        println("Julia forms:")
        for t in j_topos
            n, g = julia_gam(t)
            g_full = zeros(Int, 10, 10); g_full[1:n, 1:n] = g
            canon_sig(10, 4, g_full) == sig || continue
            println("  J id=$(t[:id]): $(t[:gam_str])")
        end
        println("qgraf forms:")
        reps = get(q_sig_to_reps, sig, [])
        for e in reps
            g = qgraf_gam(e)
            edges = String[]
            for i in 1:e[:n], j in i:e[:n]
                mu = i == j ? g[i, i] : g[i, j]
                mu > 0 && push!(edges, "($i,$j)=$mu")
            end
            println("  Q (vdeg=$(e[:vdeg])):  $(join(edges, " "))")
        end
    end
end

if length(only_j) > 0
    println("\n=== Julia topology signatures NOT in qgraf ===")
    for sig in only_j
        tids = j_sig_to_tid[sig]
        println("  Julia TOPO IDs: $tids")
        for tid in tids
            t = first(filter(x -> x[:id] == tid, j_topos))
            println("    id=$tid nem=$(t[:nem]) |G|=$(t[:G]) contrib=$(t[:contrib])")
            println("      n_self=$(t[:n_self]) n_par=$(t[:n_parallel]) n_di=$(t[:n_di])")
            println("      gam=$(t[:gam_str])")
        end
    end
end

if length(only_q) > 0
    println("\n=== qgraf topology signatures NOT in Julia ===")
    for sig in only_q
        println("  $(q_topo_sigs[sig]) emissions")
    end
end
