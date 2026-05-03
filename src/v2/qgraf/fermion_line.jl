#  Phase 18b-3: multi-vertex fermion-line traversal.
#
#  walk_fermion_lines locates each fermion path through an emitted
#  diagram. A fermion line is an alternating chain
#      bar_external → vertex₁ → propagator₁ → vertex₂ → ... → vertexₙ → plain_external
#  where each internal vertex has exactly two fermion half-edges (the
#  line enters at one and exits at the other). For tree QED 2→2 the
#  line has n=1 vertex (no internal propagator); for Compton tree it
#  has n=2 vertices joined by 1 internal fermion propagator.
#
#  The bar/plain assignment uses the externals' position metadata
#  from build_externals: the line starts at a `:left` external and
#  terminates at a `:right` external.

"""
    FermionLine

One fermion line through an emitted diagram. Holds the ordered list of
internal vertices along the line plus, for each vertex, the slot indices
where the line enters (bar side) and exits (plain side). Internal
fermion propagators between consecutive vertices are identified by
their global edge ids.

Fields:
- `vertices`: ordered internal vertices along the line (length n ≥ 1).
- `in_slots[k]`: slot at `vertices[k]` where the line enters.
- `out_slots[k]`: slot at `vertices[k]` where the line exits.
- `propagator_edge_ids[k]`: global edge id of the internal fermion
  propagator joining `vertices[k]` to `vertices[k+1]` (length n − 1).
- `bar_leg, plain_leg`: external leg indices (1..n_ext) at the two ends.
  `vmap[vertices[1], in_slots[1]] == bar_leg` and
  `vmap[vertices[end], out_slots[end]] == plain_leg`.
"""
struct FermionLine
    vertices::Vector{Int}
    in_slots::Vector{Int}
    out_slots::Vector{Int}
    propagator_edge_ids::Vector{Int}
    bar_leg::Int
    plain_leg::Int
end

"""
    walk_fermion_lines(state, labels, pmap, physical_moms, n_inco, model;
                       ps1, phys_anti=nothing) -> Vector{FermionLine}

Walk each fermion line from its `:left` external endpoint through any
internal vertices and propagators until reaching a `:right` external
endpoint. Each internal vertex on a line must have exactly two fermion
half-edges (QED/QCD/EW 3-vertices satisfy this).

`ps1` and `phys_anti` are threaded to `build_externals` so each slot's
bar/plain position reflects the PHYSICAL particle/antiparticle identity
of the ps1-permuted leg.
"""
function walk_fermion_lines(state::TopoState, labels,
                             pmap::AbstractMatrix{Symbol},
                             physical_moms::Vector{Momentum},
                             n_inco::Int,
                             model::AbstractModel;
                             ps1::AbstractVector{<:Integer}=1:Int(state.n_ext),
                             phys_anti::Union{Nothing, Vector{Bool}}=nothing)
    n_ext = Int(state.n_ext)
    externals = build_externals(state, pmap, physical_moms, n_inco, model;
                                  ps1=ps1, phys_anti=phys_anti)
    amap = compute_amap(state, labels)

    out = FermionLine[]
    for left_leg in 1:n_ext
        externals[left_leg].position == :left || continue

        # External legs have a single slot (slot 1). Step into the
        # internal vertex they attach to.
        v_curr  = Int(labels.vmap[left_leg, 1])
        in_slot = Int(labels.lmap[left_leg, 1])
        v_curr > n_ext ||
            error("walk_fermion_lines: external left leg $left_leg connects to another external (impossible for tree)")
        _field_species(model, pmap[v_curr, in_slot]) isa Fermion ||
            error("walk_fermion_lines: external left leg $left_leg connects to non-fermion slot at v=$v_curr s=$in_slot")

        vertices = Int[]
        in_slots = Int[]
        out_slots = Int[]
        prop_edges = Int[]

        while true
            out_slot = _other_fermion_slot(v_curr, in_slot, state, model, pmap)
            push!(vertices, v_curr)
            push!(in_slots, in_slot)
            push!(out_slots, out_slot)

            next_v = Int(labels.vmap[v_curr, out_slot])
            if next_v <= n_ext
                externals[next_v].position == :right ||
                    error("walk_fermion_lines: line from left leg $left_leg ends at non-:right external $next_v")
                push!(out, FermionLine(vertices, in_slots, out_slots,
                                        prop_edges, left_leg, next_v))
                break
            end
            push!(prop_edges, Int(amap[v_curr, out_slot]))
            in_slot = Int(labels.lmap[v_curr, out_slot])
            v_curr  = next_v
        end
    end
    out
end

# Unique fermion slot at v that is not `in_slot`. QED/QCD/EW 3-vertices
# carry exactly two fermion half-edges; higher-arity or 4-fermion
# vertices would need a different traversal policy.
function _other_fermion_slot(v::Int, in_slot::Int, state::TopoState,
                              model::AbstractModel,
                              pmap::AbstractMatrix{Symbol})
    vdeg_v = Int(state.vdeg[v])
    fermion_slots = Int[s for s in 1:vdeg_v
                          if _field_species(model, pmap[v, s]) isa Fermion]
    length(fermion_slots) == 2 ||
        error("walk_fermion_lines: vertex $v has $(length(fermion_slots)) fermion slots; only 2-fermion 3-vertices supported")
    in_slot in fermion_slots ||
        error("walk_fermion_lines: in_slot $in_slot at v=$v is not a fermion slot")
    fermion_slots[1] == in_slot ? fermion_slots[2] : fermion_slots[1]
end
