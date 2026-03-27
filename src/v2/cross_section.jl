# Layer 6: Cross-section evaluation.
#
# Design patterns:
#   - DiffEq Problem pattern (immutable problem → solve)
#   - Phase space as a callable (ps(cosθ) → Mandelstam variables)
#   - Clean separation: |M|² computation vs phase-space integration

# ---- Mandelstam kinematics for 2→2 massless ----
struct Mandelstam
    s::Rational{Int}
    t::Rational{Int}
    u::Rational{Int}
end

function Mandelstam(s::Number, cosθ::Number)
    s_r = Rational{Int}(s)
    t = -s_r // 2 * (1 - cosθ)
    u = -s_r // 2 * (1 + cosθ)
    Mandelstam(s_r, t, u)
end

# Build SPContext from Mandelstam (massless 2→2)
function sp_context_from_mandelstam(man::Mandelstam)
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
struct CrossSectionProblem
    model::AbstractModel
    incoming::Vector{ExternalLeg}
    outgoing::Vector{ExternalLeg}
    s::Number  # CM energy squared
end

# ---- Solve: the full pipeline ----
function solve_tree(prob::CrossSectionProblem)
    rules = feynman_rules(prob.model)
    diagrams = tree_diagrams(prob.model, prob.incoming, prob.outgoing)

    # For each diagram, build amplitude chains
    # Then spin-sum, trace, contract, evaluate
    all_chains = Tuple{DiracChain, DiracChain}[]
    for d in diagrams
        chains = build_amplitude(d, rules)
        push!(all_chains, chains)
    end

    # For single s-channel: spin sum → trace
    chain_e, chain_mu = all_chains[1]
    m_squared = spin_sum_amplitude_squared(chain_e, chain_mu)

    # Contract Lorentz indices
    contracted = contract(m_squared)

    # Expand and simplify
    expanded = expand_scalar_product(contracted)

    (amplitude_squared=expanded, diagrams=diagrams, rules=rules)
end

# Evaluate |M|² at a specific kinematic point
function evaluate_m_squared(result, man::Mandelstam)
    ctx = sp_context_from_mandelstam(man)
    evaluated = evaluate_sp(result; ctx=ctx)
    # Extract scalar value
    scalar_key = FactorKey()
    haskey(evaluated.terms, scalar_key) || return 0//1
    evaluated.terms[scalar_key]
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
