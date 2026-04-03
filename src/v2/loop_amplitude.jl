# Layer 3c: Loop amplitude construction for 1-loop box diagrams.
#
# Builds the numerator of the loop integral as DiracExpr chains.
# The loop momentum q is a regular Momentum — MomentumSum handles q+p₁ etc.
# The denominator structure (propagator momenta and masses) is returned
# separately for use by the tensor integral decomposition (tid.jl).
#
# Momentum routing convention (Denner, Eq. 4.1):
#   D₀ = q²,  D₁ = (q+K₁)²-m₁²,  D₂ = (q+K₂)²-m₂²,  D₃ = (q+K₃)²-m₃²
# where K₁, K₂, K₃ are accumulated external momenta around the loop.
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 4.4, Eq. (4.1)
# Ref: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/ElAel-MuAmu.m

"""
    BoxDenominators

Propagator structure of a box integral in PaVe convention.
Stores accumulated momenta and the Mandelstam-dependent invariant labels
(:s_inv and :t_inv or :u_inv) for evaluation at specific kinematics.
"""
struct BoxDenominators
    accum_mom_names::NTuple{3, Any}         # K₁, K₂, K₃ as Momentum/MomentumSum
    masses::NTuple{4, Rational{Int}}        # m₀², m₁², m₂², m₃²
    kinematic_inv::Symbol                   # :t or :u — which Mandelstam is p₁₃
end

"""
    build_loop_box_amplitude(ch::LoopChannel, rules, model)

Build the 1-loop box amplitude numerator. Returns:
- `chain_e::DiracExpr` — electron fermion line with loop momentum
- `chain_mu::DiracExpr` — muon fermion line with loop momentum
- `denoms::BoxDenominators` — denominator structure for TID

Ref: Denner1993, Eq. (4.1): "the amplitude is ∫ d^D q N(q)/[D₀D₁D₂D₃]"
"""
function build_loop_box_amplitude(ch::LoopChannel, rules::FeynmanRules,
                                  model::AbstractModel)
    in1, in2, out1, out2 = ch.legs
    q = Momentum(ch.loop_momentum)

    # Photon Lorentz indices (shared between the two fermion lines)
    mu = LorentzIndex(:mu_box1, DimD())
    nu = LorentzIndex(:nu_box2, DimD())

    # Electron line is the SAME for both topologies:
    # v̄(p₂) γ^ν GS(q+p₁) γ^μ u(p₁)
    # The internal electron carries momentum q+p₁.
    prop_e = MomentumSum([(1//1, q), (1//1, in1.momentum)])
    chain_e = _box_fermion_line(in1, in2, prop_e, mu, nu)

    if ch.topology == :direct_box
        # Direct: photon μ connects e⁻ to μ⁻, photon ν connects e⁺ to μ⁺.
        # Muon line: ū(k₁) γ_μ GS(q+k₁) γ_ν v(k₂)
        # Internal muon carries momentum q+k₁.
        prop_mu = MomentumSum([(1//1, q), (1//1, out1.momentum)])
        chain_mu = _box_fermion_line(out1, out2, prop_mu, mu, nu)
        denoms = BoxDenominators(
            (in1.momentum,
             MomentumSum([(1//1, in1.momentum), (1//1, in2.momentum)]),
             out1.momentum),
            (0//1, 0//1, 0//1, 0//1),  # all massless
            :t)  # p₁₃ = (p₁-k₁)² = t
    elseif ch.topology == :crossed_box
        # Crossed: photon μ connects e⁻ to μ⁺, photon ν connects e⁺ to μ⁻.
        # Muon line: ū(k₁) γ_ν GS(q+k₂) γ_μ v(k₂)  [μ↔ν swapped, k₁→k₂]
        # Internal muon carries momentum q+k₂.
        prop_mu = MomentumSum([(1//1, q), (1//1, out2.momentum)])
        chain_mu = _box_fermion_line(out1, out2, prop_mu, nu, mu)  # note: nu,mu swapped
        denoms = BoxDenominators(
            (in1.momentum,
             MomentumSum([(1//1, in1.momentum), (1//1, in2.momentum)]),
             out2.momentum),
            (0//1, 0//1, 0//1, 0//1),
            :u)  # p₁₃ = (p₁-k₂)² = u
    else
        error("Unknown box topology: $(ch.topology)")
    end
    (chain_e, chain_mu, denoms)
end

# Build one fermion line of the box.
#
# The line has 2 vertices (γ^{mu1} and γ^{mu2}) and one internal
# fermion propagator with loop momentum (prop_mom).
#
# mu1 is the photon index at the vertex closer to leg_a (leg with
# higher spinor position in the chain). mu2 is at leg_b's vertex.
#
# Ref: Denner1993, Eq. (A.1) for QED vertex structure
function _box_fermion_line(leg_a::ExternalLeg, leg_b::ExternalLeg,
                           prop_mom::MomentumSum,
                           mu_a::LorentzIndex, mu_b::LorentzIndex)
    sp_a, pos_a = _spinor_and_position(leg_a)
    sp_b, pos_b = _spinor_and_position(leg_b)
    bar_sp = pos_a == :left ? sp_a : sp_b
    plain_sp = pos_a == :left ? sp_b : sp_a

    gs_prop = DiracGamma(MomSumSlot(prop_mom))
    g_mu_a = DiracGamma(LISlot(mu_a))
    g_mu_b = DiracGamma(LISlot(mu_b))

    # Chain: bar_sp × [γ at bar's vertex] × GS(prop) × [γ at plain's vertex] × plain_sp
    if pos_a == :left
        # leg_a is bar (left), leg_b is plain (right)
        # bar vertex gets mu_a, plain vertex gets mu_b
        chain = dot(bar_sp, g_mu_a, gs_prop, g_mu_b, plain_sp)
    else
        # leg_a is plain (right), leg_b is bar (left)
        # bar vertex gets mu_b, plain vertex gets mu_a
        chain = dot(bar_sp, g_mu_b, gs_prop, g_mu_a, plain_sp)
    end
    DiracExpr(chain)
end
