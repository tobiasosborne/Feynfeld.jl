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
                           physical_moms, n_inco) -> AmplitudeBundle

Master assembler for one emission. Composes route_momenta + compute_amap
+ build_propagators + build_vertices + build_externals + walk_fermion_lines.

Convention:
- `physical_moms` holds the PHYSICAL leg momenta in slot order
  (incoming first, then outgoing — same order as `pmap[1..n_ext, 1]`).
- `n_inco` is the number of incoming legs.
- Internally builds the qgraf "all incoming" ext_moms by negating the
  outgoing-leg momenta for momentum routing.

Tree-only scope; the helper bails with a Phase-18b deferral message
when `walk_fermion_lines` rejects an internal fermion propagator.
"""
function emission_to_amplitude(state::TopoState, labels,
                                ps1::AbstractVector{<:Integer},
                                pmap::AbstractMatrix{Symbol},
                                model::AbstractModel;
                                physical_moms::Vector{Momentum},
                                n_inco::Int)
    n_ext = Int(state.n_ext)
    length(physical_moms) == n_ext ||
        error("emission_to_amplitude: physical_moms length $(length(physical_moms)) ≠ n_ext $n_ext")

    # qgraf "all incoming" momenta: keep incoming as-is, negate outgoing.
    qgraf_ext_moms = Vector{Momentum}(undef, n_ext)
    for i in 1:n_ext
        qgraf_ext_moms[i] = physical_moms[i]
    end

    edge_mom    = route_momenta(state, labels, qgraf_ext_moms)
    propagators = build_propagators(state, labels, pmap, edge_mom, model)
    vertices    = build_vertices(state, labels, pmap, edge_mom, model)
    externals   = build_externals(state, pmap, physical_moms, n_inco, model)
    lines       = walk_fermion_lines(state, labels, pmap, physical_moms,
                                       n_inco, model)

    # Per-line chain: bar_spinor × vertex_factor × plain_spinor.
    # Mirror of src/v2/amplitude.jl:61-72 _fermion_line_chain.
    line_chains = DiracExpr[]
    for line in lines
        bar_sp   = externals[line.bar_leg].spinor
        plain_sp = externals[line.plain_leg].spinor
        vtx      = vertices[line.vertex]
        chain_terms = Tuple{AlgSum, Main.FeynfeldX.DiracChain}[]
        for (coeff, chain) in vtx.terms
            full = Main.FeynfeldX.dot(bar_sp, chain.elements..., plain_sp)
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
    sym_factor   = Rational{Int}(Main.FeynfeldX.QgrafPort.compute_local_sym_factor(
        state, labels, pmap, _conjugate_dict(model)))

    AmplitudeBundle(line_chains, amplitude, denoms,
                    fermion_sign, 1 // sym_factor, alg(1))
end

# Build conjugate dict on demand from the model (mirrors
# diagram_gen.jl::_expand_model_for_diagen).
function _conjugate_dict(model::AbstractModel)
    d = Dict{Symbol, Symbol}()
    for f in Main.FeynfeldX.model_fields(model)
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
    for f in Main.FeynfeldX.model_fields(model)
        v = f isa Main.FeynfeldX.Field{Fermion} ? 1 : 0
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
