#  Phase 18a-7: master assembler emission_to_amplitude.
#
#  Composes Phase 18a-1..6 into an AmplitudeBundle for one emission.
#  Tree-only scope; multi-vertex fermion lines (Compton tree, fermion
#  loops) and explicit boson polarisation are deferred to Phase 18b.
#
#  qgraf "all incoming" convention is used internally for momentum
#  routing; physical momenta are passed separately for spinor lookups.

"""
    AmplitudeBundle

Per-emission amplitude bundle.

Fields:
- `line_chains`: Vector{DiracExpr} — one per fermion line. For the
  existing `spin_sum_amplitude_squared(de1, de2)` interface (Layer
  4b) we keep them separate; tree QED 2→2 has 2 lines, φ³ has 0.
- `amplitude`: DiracExpr — convenience: product of `line_chains`
  (or `DiracExpr(alg(1))` when there are no lines).
- `denoms`: Vector{AlgSum} — one (p² − m²) factor per internal
  propagator, to be inverted by the cross-section evaluator.
- `fermion_sign`: ±1 from qdis_fermion_sign.
- `sym_factor`: 1/S_local from compute_local_sym_factor (Rational).
- `coupling`: AlgSum placeholder for the coupling product (e.g. e²);
  Phase 18a fills with alg(1), full coupling assignment is 18b.
"""
struct AmplitudeBundle
    line_chains::Vector{DiracExpr}
    amplitude::DiracExpr
    denoms::Vector{AlgSum}
    fermion_sign::Int
    sym_factor::Rational{Int}
    coupling::AlgSum
end

"""
    emission_to_amplitude(state, labels, ps1, pmap, model;
                           physical_moms, n_inco, phys_anti=nothing)
        -> AmplitudeBundle

Master assembler for one emission. Composes route_momenta + compute_amap
+ build_propagators + build_vertices + build_externals + walk_fermion_lines.

Convention:
- `physical_moms` holds the PHYSICAL leg momenta in PHYSICAL-leg order
  (incoming 1..n_inco, outgoing n_inco+1..n_ext).
- `n_inco` is the number of incoming legs.
- `ps1[i]` is the physical-leg index at qgen slot i (forward permutation).
  Slot i sees the physical momentum `physical_moms[ps1[i]]` and inherits
  physical incoming/outgoing status from `ps1[i] ≤ n_inco`.
- For `route_momenta` we use the qgraf "all incoming" convention: slots
  whose physical leg is outgoing (`ps1[i] > n_inco`) contribute with a
  −1 sign via `ext_signs`. Mirrors qgraf's `qflow` output-time sign flip
  (f08:6961-6964). Spinors receive unnegated physical momenta.
- `phys_anti[i]` (optional) is the PHYSICAL antiparticle flag of leg i
  (in physical-leg order). When supplied, threaded to `build_externals`
  / `walk_fermion_lines` and used by `_spinor_dispatch` instead of the
  qgraf-pmap-label-derived flag. Required for Bhabha-class processes
  where a slot's pmap label (qgraf all-incoming convention) disagrees
  with the leg's physical particle/antiparticle identity. Defaults to
  `nothing` (back-compat: dispatch falls back to label-derived).

Tree-only scope; the helper bails with a Phase-18b deferral message
when `walk_fermion_lines` rejects an internal fermion propagator.
"""
function emission_to_amplitude(state::TopoState, labels,
                                ps1::AbstractVector{<:Integer},
                                pmap::AbstractMatrix{Symbol},
                                model::AbstractModel;
                                physical_moms::Vector{Momentum},
                                n_inco::Int,
                                phys_anti::Union{Nothing, Vector{Bool}}=nothing)
    n_ext = Int(state.n_ext)
    length(physical_moms) == n_ext ||
        error("emission_to_amplitude: physical_moms length $(length(physical_moms)) ≠ n_ext $n_ext")
    length(ps1) == n_ext ||
        error("emission_to_amplitude: ps1 length $(length(ps1)) ≠ n_ext $n_ext")
    phys_anti === nothing || length(phys_anti) == n_ext ||
        error("emission_to_amplitude: phys_anti length $(length(phys_anti)) ≠ n_ext $n_ext")

    # Slot i gets physical leg ps1[i]; outgoing slots flip sign for
    # qgraf "all incoming" routing convention.
    qgraf_ext_moms = Vector{Momentum}(undef, n_ext)
    ext_signs      = Vector{Int}(undef, n_ext)
    for i in 1:n_ext
        phys_idx          = Int(ps1[i])
        qgraf_ext_moms[i] = physical_moms[phys_idx]
        ext_signs[i]      = phys_idx > n_inco ? -1 : 1
    end

    edge_mom    = route_momenta(state, labels, qgraf_ext_moms; ext_signs=ext_signs)
    propagators = build_propagators(state, labels, pmap, edge_mom, model)
    vertices    = build_vertices(state, labels, pmap, edge_mom, model)
    externals   = build_externals(state, pmap, physical_moms, n_inco, model;
                                    ps1=ps1, phys_anti=phys_anti)
    lines       = walk_fermion_lines(state, labels, pmap, physical_moms,
                                       n_inco, model;
                                       ps1=ps1, phys_anti=phys_anti)

    # Per-line chain: bar_spinor × vertex_factor × plain_spinor.
    # Mirror of src/v2/amplitude.jl:61-72 _fermion_line_chain.
    line_chains = DiracExpr[]
    for line in lines
        bar_sp   = externals[line.bar_leg].spinor
        plain_sp = externals[line.plain_leg].spinor
        vtx      = vertices[line.vertex]
        chain_terms = Tuple{AlgSum, DiracChain}[]
        for (coeff, chain) in vtx.terms
            full = dot(bar_sp, chain.elements..., plain_sp)
            push!(chain_terms, (coeff, full))
        end
        push!(line_chains, DiracExpr(chain_terms))
    end

    # Master amplitude: product of all line chains. For φ³ (no lines)
    # the amplitude is just DiracExpr(alg(1)).
    amplitude = isempty(line_chains) ? DiracExpr(alg(1)) :
                foldl(*, line_chains)

    denoms = AlgSum[p.denom for p in propagators]

    fermion_sign = _emission_fermion_sign(state, labels, pmap, ps1, n_inco, model)
    sym_factor   = Rational{Int}(compute_local_sym_factor(
        state, labels, pmap, _conjugate_dict(model)))

    AmplitudeBundle(line_chains, amplitude, denoms,
                    fermion_sign, 1 // sym_factor, alg(1))
end

# Build conjugate dict on demand from the model (mirrors
# diagram_gen.jl::_expand_model_for_diagen).
function _conjugate_dict(model::AbstractModel)
    d = Dict{Symbol, Symbol}()
    for f in model_fields(model)
        if f.self_conjugate
            d[f.name] = f.name
        else
            anti = Symbol(f.name, :_bar)
            d[f.name] = anti
            d[anti]   = f.name
        end
    end
    d
end

function _antiq_dict(model::AbstractModel)
    d = Dict{Symbol, Int}()
    for f in model_fields(model)
        v = f isa Field{Fermion} ? 1 : 0
        d[f.name] = v
        f.self_conjugate || (d[Symbol(f.name, :_bar)] = v)
    end
    d
end

function _emission_fermion_sign(state, labels, pmap, ps1, n_inco, model)
    amap = compute_amap(state, labels)
    qdis_fermion_sign(state, labels, pmap, ps1, n_inco,
                      _antiq_dict(model), _conjugate_dict(model), amap)
end
