# Submodule root for Strategy C qgraf port.
# Isolates the new types from the legacy FeynmanTopology / DegreePartition.
# Ref: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08
# Ref: refs/papers/Nogueira1993_JCompPhys105_279.pdf
module QgrafPort

include("types.jl")
include("canonical.jl")
include("topology.jl")
include("qgen.jl")

export Partition, EquivClass, FilterSet, TopoState,
       rho_k, n_internal, n_vertices, n_edges, no_filters,
       MAX_V,
       compute_equiv_classes!, next_class_perm!, is_canonical_full!,
       is_canonical_feynman,
       step_b_enumerate!,
       step_c_enumerate!,
       qg21_enumerate!,
       qg10_enumerate!,
       _is_connected_internal,
       build_dpntro,
       compute_qg10_labels

end
