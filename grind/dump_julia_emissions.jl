#  Dump Julia qgen emissions for QED ee→μμ 1L for direct comparison with qgraf.
#  Mirrors the GRIND format: per-emission vlis, vfo (rule index), pmap.

include("../src/v2/FeynfeldX.jl")
using .FeynfeldX
using .FeynfeldX.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              compute_qg10_labels, build_dpntro,
                              qgen_enumerate_assignments,
                              enumerate_topology_automorphisms,
                              emission_stabilizer

m = qed_model()
in_f  = [:e, :e]
out_f = [:mu, :mu]
loops = 1
n_ext = length(in_f) + length(out_f)
ext_raw = vcat(in_f, out_f)
exp_m   = FeynfeldX._expand_model_for_diagen(m)
ext_exp = FeynfeldX._expand_external_fields(ext_raw, exp_m)
dpntro  = build_dpntro(exp_m.vertex_rules)
rules   = FeynfeldX.feynman_rules(m)
vd      = Set(length(k) for k in keys(rules.vertices))

# Dump dpntro for comparison with qgraf's rotvpo
println(stderr, "JULIA DPNTRO --- field-id-free dump (we use Symbol names) ---")
for (deg, rules_for_deg) in sort(collect(dpntro))
    println(stderr, "JULIA DPNTRO deg=$deg n_rules=$(length(rules_for_deg))")
    for (i, r) in enumerate(rules_for_deg)
        println(stderr, "JULIA DPNTRO   rule deg=$deg #$i fields=$(r)")
    end
end
println(stderr, "JULIA EXT_EXP=$ext_exp")

# Iterate topologies + ext perms + qgen emissions
n_topos = Ref(0)
n_ps1 = Ref(0)
n_emissions = Ref(0)
total_burnside = Ref(Rational{Int}(0))
emit_id = Ref(0)

for dp in FeynfeldX._degree_partitions(n_ext, loops, vd)
    degs = sort([d for (d, c) in dp.counts if c > 0])
    isempty(degs) && continue
    cv = Int8[get(dp.counts, d, 0) for d in degs[1]:degs[end]]
    qp = Partition(Int8(n_ext), cv, Int8(degs[1]), Int8(loops))
    s  = TopoState(qp)
    qg21_enumerate!(s) do state
        n_topos[] += 1
        topo_id = n_topos[]
        labels = compute_qg10_labels(state)
        autos  = enumerate_topology_automorphisms(state)
        g_size = length(autos)
        # Print topology details
        println(stderr, "JULIA ===== TOPO #$topo_id =====")
        println(stderr, "JULIA TOPO n=$(state.n) n_ext=$(state.n_ext) rhop1=$(state.rhop1)")
        println(stderr, "JULIA TOPO vdeg=$(state.vdeg[1:Int(state.n)])")
        # xg upper triangle
        n = Int(state.n)
        xg_str = String[]
        for i in 1:n, j in i:n
            g = state.xg[i, j]
            if g != 0
                push!(xg_str, "($i,$j)=$g")
            end
        end
        println(stderr, "JULIA TOPO xg=", join(xg_str, " "))
        println(stderr, "JULIA TOPO vlis=$(labels.vlis[1:n])")
        println(stderr, "JULIA TOPO vmap rows: ")
        for v in 1:n
            row = labels.vmap[v, 1:Int(state.vdeg[v])]
            println(stderr, "  v$v vmap=$row lmap=$(labels.lmap[v, 1:Int(state.vdeg[v])]) rdeg=$(labels.rdeg[v]) sdeg=$(labels.sdeg[v])")
        end
        println(stderr, "JULIA TOPO #autos=$g_size")
        ps1 = collect(1:n_ext)
        while true
            n_ps1[] += 1
            ext_perm = Symbol[ext_exp[ps1[i]] for i in 1:n_ext]
            qgen_enumerate_assignments(state, labels, ext_perm, dpntro,
                                        exp_m.conjugate) do st, pmap
                n_emissions[] += 1
                emit_id[] += 1
                n_stab = emission_stabilizer(st, labels, autos, ps1, pmap)
                total_burnside[] += n_stab // g_size
                println(stderr, "JULIA ===== EMIT #$(emit_id[]) (topo #$topo_id) =====")
                println(stderr, "JULIA EMIT vlis=$(labels.vlis[1:Int(st.n)])")
                println(stderr, "JULIA EMIT ps1=$(ps1)")
                println(stderr, "JULIA EMIT n_stab=$n_stab g_size=$g_size  contrib=$(n_stab//g_size)")
                for v in 1:Int(st.n)
                    fields = pmap[v, 1:Int(st.vdeg[v])]
                    println(stderr, "JULIA EMIT pmap[$v,1..vdeg]=$fields")
                end
            end
            j = n_ext - 1
            while j >= 1 && ps1[j] >= ps1[j + 1]; j -= 1; end
            j == 0 && break
            k = n_ext
            while ps1[k] <= ps1[j]; k -= 1; end
            ps1[j], ps1[k] = ps1[k], ps1[j]
            reverse!(view(ps1, (j + 1):n_ext))
        end
    end
end

println(stderr, "JULIA SUMMARY n_topos=$(n_topos[]) n_ps1_visits=$(n_ps1[]) n_emissions=$(n_emissions[])")
println(stderr, "JULIA SUMMARY burnside_sum=$(total_burnside[]) (≈$(Float64(total_burnside[])))")
