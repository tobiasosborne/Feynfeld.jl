using Feynfeld
using Feynfeld.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              compute_qg10_labels, build_dpntro,
                              qgen_count_assignments

m = qed_model()
in_f = [:e, :e]; out_f = [:mu, :mu]; loops = 1
n_ext = length(in_f) + length(out_f)
exp = Feynfeld._expand_model_for_diagen(m)
ext_exp = Feynfeld._expand_external_fields(vcat(in_f, out_f), exp)
dpntro = build_dpntro(exp.vertex_rules)

n_topo = 0
total = 0
for dp in Feynfeld._degree_partitions(n_ext, loops,
            Set(length(k) for k in keys(Feynfeld.feynman_rules(m).vertices)))
    degs = sort([d for (d, c) in dp.counts if c > 0])
    isempty(degs) && continue
    cv = Int8[get(dp.counts, d, 0) for d in degs[1]:degs[end]]
    qp = Partition(Int8(n_ext), cv, Int8(degs[1]), Int8(loops))
    s = TopoState(qp)
    qg21_enumerate!(s) do state
        global n_topo += 1
    end
end
println("Total topologies: ", n_topo)
