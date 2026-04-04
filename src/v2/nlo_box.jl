# NLO box evaluation: assemble box contributions for 2→2 processes.
#
# Runs the full pipeline for each box channel:
#   box_channels → build_loop_box_amplitude → spin_sum_tree_loop_interference
#   → contract → expand_scalar_product → evaluate_box_integral
#
# Convention: returns the loop integral in COLLIER normalization.
# The tree amplitude structure (no coupling/propagator) is embedded in the trace.
#
# Physical Born-virtual interference (spin-averaged):
#   (1/4) Σ_spins 2Re(M_tree* × M_box) = -e⁶/(32π²s) × Im(result)
# where e² = 4πα.
#
# Derivation: the overall coupling from tree* × box × loop measure is
# purely imaginary (the i from the tree propagator survives after all
# vertex/propagator phases cancel), giving 2Re(i·A·I) = -2A·Im(I).
# Verified empirically: Im(I) is forward-backward symmetric and gives
# box/Born ~ O(α/π) as expected.  Re(I) breaks F-B symmetry and gives O(1).
#
# Ref: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/ElAel-MuAmu.m
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 4

"""
    evaluate_single_box_channel(tree_amp, ch, rules, model, sp_vals, man; mu2)

Evaluate one box channel: build amplitude → interference → TID.
Returns ComplexF64 in COLLIER normalization.
"""
function evaluate_single_box_channel(
    tree_amp::NTuple{2, DiracExpr},
    ch::LoopChannel,
    rules::FeynmanRules,
    model::AbstractModel,
    sp_vals::Dict{Tuple{Symbol,Symbol}, Float64},
    man::Mandelstam{Float64};
    mu2::Float64 = 1.0
)
    loop_e, loop_mu, denoms = build_loop_box_amplitude(ch, rules, model)
    interference = spin_sum_tree_loop_interference(tree_amp, (loop_e, loop_mu))
    contracted = contract(interference)
    expanded = expand_scalar_product(contracted)
    man_inv = if ch.topology == :direct_box
        man.t
    elseif ch.topology == :crossed_box
        man.u
    else
        error("Unknown box topology: $(ch.topology)")
    end
    evaluate_box_integral(expanded, sp_vals, denoms, man.s, man_inv; mu2=mu2)
end

"""
    evaluate_box_channels(model, rules, incoming, outgoing, man; mu2=1.0)

Sum the box loop integral over all channels for a 2→2 process.

Returns the total in COLLIER normalization (no coupling, no propagator, no spin average).
To obtain the physical spin-averaged Born-virtual interference:

    (1/4) Σ_spins 2Re(M_tree* × M_box) = -e⁶/(32π² s) × Im(result)
"""
function evaluate_box_channels(
    model::AbstractModel,
    rules::FeynmanRules,
    incoming::Vector{ExternalLeg},
    outgoing::Vector{ExternalLeg},
    man::Mandelstam{Float64};
    mu2::Float64 = 1.0
)
    sp_vals = sp_values_2to2(man)

    # Tree amplitude (first channel, s-channel for ee→μμ)
    tree_chs = tree_channels(model, rules, incoming, outgoing)
    tree_ch = first(tree_chs)
    tree_amp = build_amplitude(tree_ch, rules, model)

    # Sum over box channels
    loop_chs = box_channels(model, rules, incoming, outgoing)
    total = ComplexF64(0)
    for ch in loop_chs
        total += evaluate_single_box_channel(
            tree_amp, ch, rules, model, sp_vals, man; mu2=mu2)
    end
    total
end

"""
    born_virtual_box(model, rules, incoming, outgoing, man; alpha, mu2)

Spin-averaged box contribution to Born-virtual interference:

    (1/4) Σ_spins 2Re(M_tree* × M_box)

with coupling e² = 4πα and tree propagator 1/s.
"""
function born_virtual_box(
    model::AbstractModel,
    rules::FeynmanRules,
    incoming::Vector{ExternalLeg},
    outgoing::Vector{ExternalLeg},
    man::Mandelstam{Float64};
    alpha::Float64 = 1.0 / 137.036,
    mu2::Float64 = 1.0
)
    I_box = evaluate_box_channels(model, rules, incoming, outgoing, man; mu2=mu2)
    e_sq = 4π * alpha
    -e_sq^3 / (32π^2 * man.s) * imag(I_box)
end
