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
# Model -> Rules -> Channels -> Amplitude -> Algebra -> |M|^2
function solve_tree(prob::CrossSectionProblem)
    rules = feynman_rules(prob.model)
    channels = tree_channels(prob.model, rules, prob.incoming, prob.outgoing)
    isempty(channels) && error("No valid tree-level channels for this process")

    # Build amplitude for each channel
    all_chains = Tuple{DiracChain, DiracChain}[]
    for ch in channels
        chains = build_amplitude(ch, rules, prob.model)
        push!(all_chains, chains)
    end

    # Single channel: spin sum -> trace -> contract -> expand
    # TODO: multi-channel interference (Phase B)
    chain_L, chain_R = all_chains[1]
    m_squared = spin_sum_amplitude_squared(chain_L, chain_R)
    contracted = contract(m_squared)
    expanded = expand_scalar_product(contracted)

    (amplitude_squared=expanded, channels=channels, rules=rules)
end

# Evaluate |M|² at a specific kinematic point
function evaluate_m_squared(result, man::Mandelstam)
    ctx = sp_context_from_mandelstam(man)
    evaluated = evaluate_sp(result; ctx=ctx)
    # Extract scalar value
    scalar_key = FactorKey()
    haskey(evaluated.terms, scalar_key) || return 0//1
    coeff = evaluated.terms[scalar_key]
    coeff isa DimPoly ? evaluate_dim(coeff) : coeff
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

# ---- Compton scattering cross-section ----

# Total cross-section for Compton scattering in the massless limit.
# Ref: P&S §5.5, below Eq. (5.87). In the massless limit with s+t+u=0:
# |M̄|² = -2e⁴(u/s + s/u)
# Ref: refs/FeynCalc/.../ElGa-ElGa.m, line 96: "ampSquared[0]/.SMP["m_e"]->0"
#
# dσ/dt = |M̄|² / (16π s²)  for 2→2 massless (all identical stats: 1/(16πs²))
# Ref: P&S Eq. (4.85): dσ = |M|²/(2s) × dΦ₂, with dΦ₂ = dt/(8πs) for massless 2→2
#
# σ_total = ∫ dt × |M̄|²/(16πs²)  with t ∈ [-s, 0] (massless: u = -s-t)
function sigma_total_compton_massless(s::Float64; alpha::Float64 = 1 / 137.036)
    # |M̄|² = -2e⁴(u/s + s/u) = 2α²(4π)²(-u/s - s/u)  [e² = 4πα]
    # After substituting u = -s-t and integrating dt from -s to 0:
    # σ = (2πα²/s) × ∫₋₁⁰ dx [-(1+x)/1 - 1/(1+x)]   where x = t/s ∈ [-1,0], u/s = -1-x
    # But u/s = (-s-t)/s = -1-t/s diverges at t=0 (forward singularity).
    # This is the well-known IR divergence of massless Compton.
    # The massive case (Klein-Nishina) is IR-finite. Massless total σ is divergent.
    # For a finite result, one needs a t-cut or the full massive formula.
    error("Massless Compton total cross-section is IR divergent (forward singularity at t=0). Use massive Klein-Nishina formula.")
end
