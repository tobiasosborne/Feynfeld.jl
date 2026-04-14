#  Dump Julia qgen emissions for phi3 2L φφ→φφ for per-topology comparison
#  against qgraf-grind (grind_phi3_2L.txt).
#
#  Prints per-topology Burnside contribution so we can compare totals.

include("../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              compute_qg10_labels, build_dpntro,
                              qgen_enumerate_assignments,
                              enumerate_topology_automorphisms,
                              emission_stabilizer

m = phi3_model()
in_f  = [:phi, :phi]
out_f = [:phi, :phi]
loops = 2
n_ext = length(in_f) + length(out_f)
ext_raw = vcat(in_f, out_f)
exp_m   = FeynfeldX._expand_model_for_diagen(m)
ext_exp = FeynfeldX._expand_external_fields(ext_raw, exp_m)
dpntro  = build_dpntro(exp_m.vertex_rules)
rules   = FeynfeldX.feynman_rules(m)
vd      = Set(length(k) for k in keys(rules.vertices))

# One Julia-side emission record for later grouping.
struct EmRec
    topo_id::Int
    ps1::Vector{Int}
    gam_sig::Any      # the xg-upper-triangle as a topology signature
    g_size::Int
    n_stab::Int
    contrib::Rational{Int}
end

ems = EmRec[]
topos = Dict{Any, Int}()     # signature → topo_id (stable across emissions)
topo_stats = Dict{Int, NamedTuple}()   # topo_id → (n_emissions, sum_contrib, g_size, n_self, n_parallel, gam_str)

function _gam_sig(state)
    n = Int(state.n)
    pairs = Tuple{Int,Int,Int}[]
    for i in 1:n, j in i:n
        g = state.xg[i, j]
        if g != 0; push!(pairs, (i, j, Int(g))); end
    end
    tuple(pairs...)
end

function _gam_str(state)
    n = Int(state.n)
    out = String[]
    for i in 1:n, j in i:n
        g = state.xg[i, j]
        if g != 0; push!(out, "($i,$j)=$g"); end
    end
    join(out, " ")
end

function _topo_features(state)
    n = Int(state.n); ne = Int(state.n_ext)
    n_self = 0; n_parallel = 0; n_di = 0; n_tri = 0
    for i in (ne+1):n
        n_self += Int(state.xg[i, i])
    end
    for i in 1:n, j in (i+1):n
        g = Int(state.xg[i, j])
        if g >= 2
            n_parallel += 1
            n_di  += (g == 2)
            n_tri += (g == 3)
        end
    end
    (; n_self, n_parallel, n_di, n_tri)
end

topo_counter = Ref(0)

for dp in FeynfeldX._degree_partitions(n_ext, loops, vd)
    degs = sort([d for (d, c) in dp.counts if c > 0])
    isempty(degs) && continue
    cv = Int8[get(dp.counts, d, 0) for d in degs[1]:degs[end]]
    qp = Partition(Int8(n_ext), cv, Int8(degs[1]), Int8(loops))
    s  = TopoState(qp)
    qg21_enumerate!(s) do state
        sig    = _gam_sig(state)
        topo_id = get!(topos, sig) do
            topo_counter[] += 1
            topo_counter[]
        end
        labels = compute_qg10_labels(state)
        autos  = enumerate_topology_automorphisms(state)
        g_size = length(autos)
        features = _topo_features(state)
        ps1 = collect(1:n_ext)
        sum_contrib = Rational{Int}(0)
        n_em = 0
        while true
            ext_perm = Symbol[ext_exp[ps1[i]] for i in 1:n_ext]
            qgen_enumerate_assignments(state, labels, ext_perm, dpntro,
                                        exp_m.conjugate) do st, pmap
                n_em += 1
                n_stab = emission_stabilizer(st, labels, autos, ps1, pmap)
                contrib = n_stab // g_size
                sum_contrib += contrib
                push!(ems, EmRec(topo_id, copy(ps1), sig, g_size, n_stab, contrib))
            end
            j = n_ext - 1
            while j >= 1 && ps1[j] >= ps1[j + 1]; j -= 1; end
            j == 0 && break
            k = n_ext
            while ps1[k] <= ps1[j]; k -= 1; end
            ps1[j], ps1[k] = ps1[k], ps1[j]
            reverse!(view(ps1, (j + 1):n_ext))
        end
        # update per-topo stats (may update multiple times if qg21 revisits
        # the same canonical sig — but it shouldn't)
        topo_stats[topo_id] = (;
            n_emissions = n_em,
            sum_contrib,
            g_size,
            features...,
            gam_str = _gam_str(state),
        )
    end
end

# Print summary
ordered = sort(collect(topo_stats); by = p -> (-p.second.n_emissions, p.first))
let total_em = 0, total_contrib = Rational{Int}(0)
    for (tid, s) in ordered
        println("JULIA TOPO #$tid  #emissions=$(s.n_emissions)  |G|=$(s.g_size)  burnside_contrib=$(s.sum_contrib)")
        println("  features=(n_self=$(s.n_self), n_parallel=$(s.n_parallel), n_di=$(s.n_di), n_tri=$(s.n_tri))")
        println("  gam=$(s.gam_str)")
        total_em += s.n_emissions
        total_contrib += s.sum_contrib
    end
    println()
    println("JULIA SUMMARY: total_topos=$(length(topo_stats)) total_emissions=$(total_em) total_burnside=$(total_contrib)")
end
