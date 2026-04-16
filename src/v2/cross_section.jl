# Layer 6: Cross-section evaluation.
#
# Design patterns:
#   - DiffEq Problem pattern (immutable problem → solve)
#   - Phase space as a callable (ps(cosθ) → Mandelstam variables)
#   - Clean separation: |M|² computation vs phase-space integration

# ---- Mandelstam kinematics for 2→2 ----
struct Mandelstam{T<:Real}
    s::T
    t::T
    u::T
end

# Parametric constructor: Rational path stays exact, Float64 path avoids overflow.
function Mandelstam(s::T, cosθ::T) where {T<:Real}
    half_s = s / 2
    Mandelstam(s, -half_s * (1 - cosθ), -half_s * (1 + cosθ))
end
# Mixed-type convenience: promote all args to a common type.
Mandelstam(s::Number, cosθ::Number) = Mandelstam(promote(s, cosθ)...)
Mandelstam(s::Number, t::Number, u::Number) = Mandelstam(promote(s, t, u)...)

# Build SPContext from Mandelstam (massless 2→2, exact Rational path)
function sp_context_from_mandelstam(man::Mandelstam{Rational{Int}})
    sp_context(
        (:p1, :p1) => 0//1, (:p2, :p2) => 0//1,
        (:k1, :k1) => 0//1, (:k2, :k2) => 0//1,
        (:p1, :p2) => man.s // 2,
        (:k1, :k2) => man.s // 2,
        (:p1, :k1) => -man.t // 2,
        (:p2, :k2) => -man.t // 2,
        (:p1, :k2) => -man.u // 2,
        (:p2, :k1) => -man.u // 2,
    )
end

# ---- Cross-section problem (DiffEq pattern) ----
struct CrossSectionProblem{T<:Number}
    model::AbstractModel
    incoming::Vector{ExternalLeg}
    outgoing::Vector{ExternalLeg}
    s::T  # CM energy squared
end

# ---- Solve: the full pipeline ----
# Model -> Rules -> Channels -> Amplitude -> Algebra -> |M|^2
function solve_tree(prob::CrossSectionProblem)
    rules = feynman_rules(prob.model)
    channels = tree_channels(prob.model, rules, prob.incoming, prob.outgoing)
    isempty(channels) && error("No valid tree-level channels for this process")

    # Build amplitude for each channel
    all_amps = Tuple{DiracExpr, DiracExpr}[]
    for ch in channels
        amp = build_amplitude(ch, rules, prob.model)
        push!(all_amps, amp)
    end

    # Single channel: spin sum -> trace -> contract -> expand
    # TODO: multi-channel interference (Phase B)
    amp_L, amp_R = all_amps[1]
    m_squared = spin_sum_amplitude_squared(amp_L, amp_R)
    contracted = contract(m_squared)
    expanded = expand_scalar_product(contracted)

    (amplitude_squared=expanded, channels=channels, rules=rules)
end

# Evaluate |M|² at a specific kinematic point (exact Rational path)
function evaluate_m_squared(result, man::Mandelstam{Rational{Int}})
    ctx = sp_context_from_mandelstam(man)
    evaluated = evaluate_sp(result; ctx=ctx)
    # Extract scalar value
    scalar_key = FactorKey()
    haskey(evaluated.terms, scalar_key) || return 0//1
    coeff = evaluated.terms[scalar_key]
    coeff isa DimPoly ? evaluate_dim(coeff) : coeff
end

# ---- Float64 evaluation (bypasses SPContext, avoids Rational overflow) ----

# Build Float64 SP values for massive 2→2.
# Ref: p_i · p_j = (m_i² + m_j² - t_ij) / 2  for the appropriate Mandelstam invariant.
function sp_values_2to2(man::Mandelstam{Float64};
                        m1_sq=0.0, m2_sq=0.0, m3_sq=0.0, m4_sq=0.0)
    Dict{Tuple{Symbol,Symbol}, Float64}(
        _sp_key(:p1, :p1) => m1_sq,  _sp_key(:p2, :p2) => m2_sq,
        _sp_key(:k1, :k1) => m3_sq,  _sp_key(:k2, :k2) => m4_sq,
        _sp_key(:p1, :p2) => (man.s - m1_sq - m2_sq) / 2,
        _sp_key(:k1, :k2) => (man.s - m3_sq - m4_sq) / 2,
        _sp_key(:p1, :k1) => (m1_sq + m3_sq - man.t) / 2,
        _sp_key(:p1, :k2) => (m1_sq + m4_sq - man.u) / 2,
        _sp_key(:p2, :k1) => (m2_sq + m3_sq - man.u) / 2,
        _sp_key(:p2, :k2) => (m2_sq + m4_sq - man.t) / 2,
    )
end

# Evaluate an AlgSum numerically given Float64 scalar product values.
# Handles Pair{Momentum,Momentum} and Eps with all-Momentum slots via dispatch.
# See eps_evaluate.jl for _eval_factor dispatch methods.
function evaluate_numeric(result::AlgSum,
                          sp_vals::Dict{Tuple{Symbol,Symbol}, Float64})::Float64
    total = 0.0
    for (fk, c) in result.terms
        val = Float64(evaluate_dim(c))
        for f in fk.factors
            val *= _eval_factor(f, sp_vals)
        end
        total += val
    end
    total
end

# Float64 dispatch: bypasses SPContext, pure Float64 arithmetic.
function evaluate_m_squared(result::AlgSum, man::Mandelstam{Float64};
                            m1_sq=0.0, m2_sq=0.0, m3_sq=0.0, m4_sq=0.0)
    sp_vals = sp_values_2to2(man; m1_sq, m2_sq, m3_sq, m4_sq)
    evaluate_numeric(result, sp_vals)
end

# ---- Differential cross-section ----
# dσ/dΩ = |M|² / (64π²s)  for massless 2→2
function dsigma_domega(m_squared::Number, s::Number)
    Float64(m_squared) / (64 * π^2 * Float64(s))
end

# Total cross-section for e+e- → μ+μ- (tree, massless, analytical)
# σ = 4πα²/(3s) — P&S Eq (5.12)
# |M|² = e⁴ × 8(t²+u²) / s² (including coupling)
# After angular integration: σ = πα²/(2s) × ∫₋₁¹ (1+cos²θ) dcosθ = 4πα²/(3s)
function sigma_total_tree_ee_mumu(s::Float64; alpha::Float64=1/137.036)
    4π * alpha^2 / (3 * s)
end

# ---- Phase 18a-8: pipeline-driven tree solver ----
#
# solve_tree_pipeline drives the qg21 → qgen emission stream, builds an
# AmplitudeBundle per emission, picks the canonical orbit
# representative, and runs the standard spin-sum/contract/expand_sp
# pipeline. Existing solve_tree is preserved unchanged for back-compat.

"""
    solve_tree_pipeline(prob::CrossSectionProblem)

Build the tree-level |M|² for the given process by routing through the
qg21 diagram-emission pipeline (Phase 18a-7's emission_to_amplitude +
Phase 18b-1's Burnside combine), then running contract / expand_sp.

Returns a NamedTuple `(amplitude_squared, n_emissions, orbit_denoms)`.
`amplitude_squared` is the trace-only |M|² — sum of per-orbit diagonal
and per-orbit-pair interference traces, Burnside-weighted, with
fermion-sign factors (see src/v2/qgraf/burnside_combine.jl).
Per-orbit propagator denominators live in `orbit_denoms`; applying
the 1/denom factors is the caller's responsibility until bead
feynfeld-rj1l lands a symbolic-inverse factor (Option B retrofit).

Phase 18b-1 scope: 2-fermion-line bundles (QED/QCD tree 2→2 via boson
exchange) + 0-line bundles (φ³). Multi-vertex fermion lines are
deferred to Phase 18b-3.
"""
function solve_tree_pipeline(prob::CrossSectionProblem)
    in_fields  = [leg.field_name for leg in prob.incoming]
    out_fields = [leg.field_name for leg in prob.outgoing]
    n_inco     = length(prob.incoming)
    physical_moms = vcat([leg.momentum for leg in prob.incoming],
                         [leg.momentum for leg in prob.outgoing])

    # Phase 18b-1: keep one canonical representative per orbit.  Retaining
    # all orbit-members would require inter-bundle momentum-label matching
    # inside spin_sum_interference (interference.jl:102 `_find_line_by_bar_mom`
    # keys on bar momentum names), which breaks for automorphic relabelings
    # within the same orbit.  With canonical reps only, weights collapse to
    # 1 per orbit and the double-sum in combine_m_squared_burnside gives
    # each (orbit_i, orbit_j) pair exactly once.
    bundles = QgrafPort.AmplitudeBundle[]
    weights = Rational{Int}[]
    QgrafPort._foreach_emission(prob.model, in_fields, out_fields; loops=0) do state, labels, ps1, pmap
        autos = QgrafPort.enumerate_topology_automorphisms(state)
        QgrafPort.is_emission_canonical(state, labels, autos, ps1, pmap) || return
        bundle = QgrafPort.emission_to_amplitude(state, labels, ps1, pmap,
                                                   prob.model;
                                                   physical_moms, n_inco)
        push!(bundles, bundle)
        push!(weights, 1//1)
    end
    isempty(bundles) && error("solve_tree_pipeline: no emissions for this process")

    m_sq       = QgrafPort.combine_m_squared_burnside(bundles, weights)
    contracted = contract(m_sq)
    expanded   = expand_scalar_product(contracted)

    (amplitude_squared=expanded,
     n_emissions=length(bundles),
     orbit_denoms=[b.denoms for b in bundles])
end
