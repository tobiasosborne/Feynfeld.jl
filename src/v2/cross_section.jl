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
