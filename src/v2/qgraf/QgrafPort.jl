# Submodule root for Strategy C qgraf port.
# Isolates the new types from the legacy FeynmanTopology / DegreePartition.
# Ref: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08
# Ref: refs/papers/Nogueira1993_JCompPhys105_279.pdf
module QgrafPort

import ..FeynfeldX: Momentum, MomentumSum, momentum_sum,
                     AlgSum, alg, pair, DiracExpr, DiracChain, dot,
                     AbstractModel, model_fields, get_field, species, Field,
                     Boson, Fermion, Scalar

include("types.jl")
include("canonical.jl")
include("topology.jl")
include("qgen.jl")
include("filters.jl")
include("spanning_tree.jl")
include("momentum.jl")
include("halfedge.jl")
include("propagator_assemble.jl")
include("vertex_assemble.jl")
include("fermion_line.jl")
include("emission_amplitude.jl")
include("audition.jl")

export Partition, EquivClass, FilterSet, TopoState,
       rho_k, n_internal, n_vertices, n_edges, no_filters,
       MAX_V,
       compute_equiv_classes!, next_class_perm!, is_canonical_full!,
       is_canonical_qgraf!, is_canonical_feynman, enumerate_topology_automorphisms,
       step_b_enumerate!,
       step_c_enumerate!,
       qg21_enumerate!,
       qg10_enumerate!,
       _is_connected_internal,
       build_dpntro,
       compute_qg10_labels,
       qgen_count_assignments,
       qdis_fermion_sign,
       compute_local_sym_factor,
       has_no_selfloop, has_no_diloop, has_no_parallel,
       is_one_pi, has_no_sbridge, has_no_tadpole, has_no_onshell,
       has_no_snail, is_one_vi,
       build_spanning_tree, count_chords,
       InternalEdge, EdgeMomenta, route_momenta,
       compute_amap,
       Propagator, build_propagators,
       build_vertices,
       ExternalFactor, build_externals,
       FermionLine, walk_fermion_lines,
       AmplitudeBundle, emission_to_amplitude,
       _foreach_emission,
       qgen_enumerate_assignments,
       count_dedup_burnside, count_dedup_canonical, count_dedup_prefilter,
       count_diagrams_qg21,
       is_emission_canonical, emission_stabilizer,
       is_pmap_canonical, pmap_stabilizer

end
