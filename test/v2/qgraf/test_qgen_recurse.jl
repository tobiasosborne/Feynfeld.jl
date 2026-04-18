#  Phase 12c: qgen recursive backtracker (multiset rule matching).
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13880-13987.
#
#  This minimal port returns RAW counts (sum over (perm, assignment) tuples).
#  It over-counts by the topology automorphism size — the symmetry-factor
#  1/S weighting that gives the true distinct-diagram count is Phase 15.

using Test
using Feynfeld
using Feynfeld: _expand_model_for_diagen, _expand_external_fields
using Feynfeld.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              compute_qg10_labels, build_dpntro,
                              qgen_count_assignments

@testset "Phase 12c: qgen recursive backtracker (raw counts)" begin

    "Driver: total field-assignment count summed over (topology, ext-perm)."
    function _qgen_total(model, in_fields::Vector{Symbol},
                          out_fields::Vector{Symbol}; loops::Int)
        n_ext = length(in_fields) + length(out_fields)
        ext_raw = vcat(in_fields, out_fields)
        exp = _expand_model_for_diagen(model)
        ext_expanded = _expand_external_fields(ext_raw, exp)
        dpntro = build_dpntro(exp.vertex_rules)

        # For Phase 12c, exhaust qg21 topologies; for each, exhaust qg10
        # ext-leg perms; for each perm, run qgen with ps1-permuted externals.
        # We bypass the qg10_enumerate!+qg21 stack here for clarity.
        rules = feynman_rules(model)
        vd = Set(length(k) for k in keys(rules.vertices))
        total = 0
        # Use the legacy partition enumerator (Phase 10 not yet ported).
        for dp in Feynfeld._degree_partitions(n_ext, loops, vd)
            degs = sort([d for (d, c) in dp.counts if c > 0])
            mrho = isempty(degs) ? 0 : degs[1]
            counts_vec = isempty(degs) ? Int8[] :
                         Int8[get(dp.counts, d, 0) for d in mrho:degs[end]]
            qp = Partition(Int8(n_ext), counts_vec, Int8(mrho), Int8(loops))
            s = TopoState(qp)
            qg21_enumerate!(s) do state
                labels = compute_qg10_labels(state)
                # Iterate all n_ext! perms of ps1
                ps1 = collect(1:n_ext)
                while true
                    ext_perm = Symbol[ext_expanded[ps1[i]] for i in 1:n_ext]
                    total += qgen_count_assignments(state, labels, ext_perm,
                                                     dpntro, exp.conjugate)
                    # Lex-next perm (Knuth Algorithm L on Vector{Int})
                    j = n_ext - 1
                    while j >= 1 && ps1[j] >= ps1[j + 1]
                        j -= 1
                    end
                    j == 0 && break
                    k = n_ext
                    while ps1[k] <= ps1[j]
                        k -= 1
                    end
                    ps1[j], ps1[k] = ps1[k], ps1[j]
                    reverse!(view(ps1, (j + 1):n_ext))
                end
            end
        end
        total
    end

    @testset "phi3 tree φ→φφ — raw count = 6 (over-counted by 3!)" begin
        # Golden master: 1 diagram.  qgen raw: 1 topo × 6 ext-perms × 1 assign = 6.
        # Phase 15 will divide by topology automorphism size (3! for 3 identical externals).
        m = phi3_model()
        @test _qgen_total(m, [:phi], [:phi, :phi]; loops=0) == 6
    end

    @testset "phi3 tree φφ→φφ — raw count = 24 (over-counted by 8)" begin
        # Golden master: 3 diagrams (s/t/u).  qgen raw: 1 topo × 24 perms × 1 assign = 24.
        # Phase 15 divisor: 8 = (2 internal swap) × (2 ext-pair-1) × (2 ext-pair-2).
        m = phi3_model()
        @test _qgen_total(m, [:phi, :phi], [:phi, :phi]; loops=0) == 24
    end

    @testset "QED ee→μμ tree — raw count = 8 (over-counted by 8)" begin
        # Golden master: 1 diagram.  Field assignment is restrictive: only
        # perms placing ee̅ together at one vertex and μμ̅ at the other yield
        # a valid assignment.  Hand count: 2(swap-vertex) × 2(swap ee̅) × 2(swap μμ̅) = 8.
        m = qed_model()
        @test _qgen_total(m, [:e, :e], [:mu, :mu]; loops=0) == 8
    end

end
