#  Phase 18a-6: fermion-line traversal.
#
#  walk_fermion_lines locates each fermion path through an emitted
#  diagram. Phase 18a tree-only scope: every fermion line touches
#  exactly ONE internal vertex (no internal fermion propagator), so
#  the FermionLine struct collapses to a single (vertex, bar_slot,
#  plain_slot, bar_leg, plain_leg) tuple. Multi-vertex fermion lines
#  (Compton tree, fermion loops) are deferred to Phase 18b.
#
#  The bar/plain assignment uses the externals' position metadata
#  from build_externals (Phase 18a-5): the external connected via
#  vmap to the bar slot must have :left position, the one at the
#  plain slot must have :right position.

"""
    FermionLine

One fermion line through an emitted diagram. Tree-only: one vertex per
line.

Fields:
- `vertex`: the (single) internal vertex on this line.
- `bar_slot, plain_slot`: slot indices at `vertex` where the line
  enters from the bar (`:left`) external and exits at the plain
  (`:right`) external.
- `bar_leg, plain_leg`: external leg indices (1..n_ext) at the two
  ends. `vmap[vertex, bar_slot] == bar_leg` and similarly for plain.
"""
struct FermionLine
    vertex::Int
    bar_slot::Int
    plain_slot::Int
    bar_leg::Int
    plain_leg::Int
end

"""
    walk_fermion_lines(state, labels, pmap, physical_moms, n_inco, model) -> Vector{FermionLine}

For each internal vertex, locate the (single) pair of fermion half-edges
and pair them into a FermionLine. Internal fermion propagators (where
the second endpoint is also internal) error with a Phase-18b deferral
message.
"""
function walk_fermion_lines(state::TopoState, labels,
                             pmap::AbstractMatrix{Symbol},
                             physical_moms::Vector{Momentum},
                             n_inco::Int,
                             model::AbstractModel)
    n_ext = Int(state.n_ext)
    rhop1 = Int(state.rhop1)
    n     = Int(state.n)
    externals = build_externals(state, pmap, physical_moms, n_inco, model)

    out = FermionLine[]
    for v in rhop1:n
        vdeg_v = Int(state.vdeg[v])
        fermion_slots = Int[s for s in 1:vdeg_v
                              if _field_species(model, pmap[v, s]) isa Fermion]
        isempty(fermion_slots) && continue

        # Tree-only: each fermion slot must connect to an external leg.
        for s in fermion_slots
            Int(labels.vmap[v, s]) <= n_ext ||
                error("walk_fermion_lines: vertex $v slot $s connects to internal vertex $(Int(labels.vmap[v,s])); multi-vertex fermion lines deferred to Phase 18b")
        end

        # Pair them: exactly one bar (left) + one plain (right) per vertex
        # for QED 3-vertices. Generalise via filtering.
        bar_slot   = 0
        plain_slot = 0
        for s in fermion_slots
            ext_leg  = Int(labels.vmap[v, s])
            position = externals[ext_leg].position
            if position == :left
                bar_slot == 0 ||
                    error("walk_fermion_lines: vertex $v has multiple bar-spinor legs (deferred to Phase 18b)")
                bar_slot = s
            elseif position == :right
                plain_slot == 0 ||
                    error("walk_fermion_lines: vertex $v has multiple plain-spinor legs (deferred to Phase 18b)")
                plain_slot = s
            end
        end
        bar_slot != 0 && plain_slot != 0 ||
            error("walk_fermion_lines: vertex $v has incomplete fermion line (bar=$bar_slot, plain=$plain_slot)")

        bar_leg   = Int(labels.vmap[v, bar_slot])
        plain_leg = Int(labels.vmap[v, plain_slot])
        push!(out, FermionLine(v, bar_slot, plain_slot, bar_leg, plain_leg))
    end
    out
end
