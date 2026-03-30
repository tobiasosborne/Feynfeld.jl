# Layer 3c: Amplitude building from TreeChannel.
#
# Dispatches on exchanged field species:
# - Boson exchange: two separate DiracChains (one per fermion line)
# - Fermion exchange: single DiracChain with propagator numerator inside
#
# Ref: SPIRAL_9_PLAN.md — "Two cases determined by exchanged field species"

"""
    build_amplitude(ch::TreeChannel, rules::FeynmanRules, model::AbstractModel)

Build the tree-level amplitude for a single channel. Returns:
- Boson exchange: `(chain_L, chain_R)` — two DiracChains for two fermion lines
- Fermion exchange: `(chain_mom, chain_mass, mass)` — propagator decomposition
"""
function build_amplitude(ch::TreeChannel, rules::FeynmanRules, model::AbstractModel)
    exch = get_field(model, ch.exchanged)
    if exch isa Field{Boson}
        (legL1, legL2), (legR1, legR2) = vertex_legs(ch)
        ferm_L = _both_fermion_legs(model, legL1, legL2)
        ferm_R = _both_fermion_legs(model, legR1, legR2)
        if ferm_L && ferm_R
            _build_boson_exchange(ch)
        elseif ferm_L || ferm_R
            _build_gauge_exchange(ch, model)
        else
            error("Pure boson scattering not yet implemented")
        end
    elseif exch isa Field{Fermion}
        _build_fermion_exchange(ch, model)
    else
        error("Scalar exchange not yet implemented")
    end
end

# ---- Boson exchange ----
# Two fermion lines connected by a virtual boson.
# Channel-specific index prevents collisions in multi-channel interference.
function _build_boson_exchange(ch::TreeChannel)
    (legL1, legL2), (legR1, legR2) = vertex_legs(ch)
    mu = LorentzIndex(Symbol(:mu_, ch.channel), DimD())
    chain_L = _fermion_line_chain(legL1, legL2, mu)
    chain_R = _fermion_line_chain(legR1, legR2, mu)
    (chain_L, chain_R)
end

# Build a DiracChain for a single fermion line at a vertex with Lorentz index mu.
function _fermion_line_chain(leg_a::ExternalLeg, leg_b::ExternalLeg, mu::LorentzIndex)
    sp_a, pos_a = _spinor_and_position(leg_a)
    sp_b, pos_b = _spinor_and_position(leg_b)
    gamma = DiracGamma(LISlot(mu))  # preserves DimD from mu
    if pos_a == :left
        dot(sp_a, gamma, sp_b)
    else
        dot(sp_b, gamma, sp_a)
    end
end

# ---- Fermion exchange ----
# Single fermion line with propagator numerator (p-slash + m) inside.
# Decomposes into momentum part + mass part:
#   A = bar_sp * gamma^nu * GS(q) * gamma^mu * sp   (momentum)
#     + m * bar_sp * gamma^nu * gamma^mu * sp         (mass)
#
# Returns (chain_mom, chain_mass, mass_val) where A = chain_mom + mass * chain_mass.
function _build_fermion_exchange(ch::TreeChannel, model::AbstractModel)
    (legL1, legL2), (legR1, legR2) = vertex_legs(ch)

    # Separate fermion vs boson legs at each vertex
    ferm_L, boson_L = _separate_fermion_boson(legL1, legL2, model)
    ferm_R, boson_R = _separate_fermion_boson(legR1, legR2, model)

    # Determine spinors from the two external fermion legs
    sp_fL, pos_fL = _spinor_and_position(ferm_L)
    sp_fR, pos_fR = _spinor_and_position(ferm_R)

    # Bar spinor (left of chain), plain spinor (right of chain)
    if pos_fL == :left
        bar_sp, plain_sp = sp_fL, sp_fR
        mu_bar = LorentzIndex(Symbol(:mu_, boson_L.momentum.name), DimD())
        mu_plain = LorentzIndex(Symbol(:mu_, boson_R.momentum.name), DimD())
    else
        bar_sp, plain_sp = sp_fR, sp_fL
        mu_bar = LorentzIndex(Symbol(:mu_, boson_R.momentum.name), DimD())
        mu_plain = LorentzIndex(Symbol(:mu_, boson_L.momentum.name), DimD())
    end

    # Propagator momentum
    q = propagator_momentum(ch)

    # Propagator mass (from exchanged fermion field)
    exch = get_field(model, ch.exchanged)
    mass_val = exch.mass == :zero ? 0//1 : 1//1  # symbolic mass → use 1 as placeholder

    # Build chains: bar_sp × gamma^{mu_bar} × [propagator] × gamma^{mu_plain} × plain_sp
    chain_mom = dot(bar_sp,
        DiracGamma(LISlot(mu_bar)),
        DiracGamma(MomSumSlot(q)),
        DiracGamma(LISlot(mu_plain)),
        plain_sp)

    chain_mass = dot(bar_sp,
        DiracGamma(LISlot(mu_bar)),
        DiracGamma(LISlot(mu_plain)),
        plain_sp)

    (chain_mom, chain_mass, mass_val)
end

# Separate a pair of legs into (fermion_leg, boson_leg) at a vertex.
function _separate_fermion_boson(a::ExternalLeg, b::ExternalLeg, model::AbstractModel)
    fa = get_field(model, a.field_name)
    fb = get_field(model, b.field_name)
    if fa isa Field{Fermion} && fb isa Field{Boson}
        (a, b)
    elseif fa isa Field{Boson} && fb isa Field{Fermion}
        (b, a)
    else
        error("Expected one fermion + one boson at vertex, got $(typeof(fa)) + $(typeof(fb))")
    end
end

# Propagator momentum for the channel.
function propagator_momentum(ch::TreeChannel)
    in1, in2, out1, out2 = ch.legs
    if ch.channel == :s
        MomentumSum([(1//1, in1.momentum), (1//1, in2.momentum)])
    elseif ch.channel == :t
        MomentumSum([(1//1, in1.momentum), (-1//1, out1.momentum)])
    else  # :u
        MomentumSum([(1//1, in1.momentum), (-1//1, out2.momentum)])
    end
end

# ---- Gauge boson exchange (fermion pair ↔ triple gauge vertex) ----
# qq̄→gg s-channel: one vertex is qqg (fermion line), the other is ggg.
# Returns (chain::DiracChain, vertex::AlgSum):
#   chain = ū(p) γ^ρ v(p')  (fermion line with propagator index ρ)
#   vertex = V_{ρ,μ₁,μ₂}(k₁, k₂, -q)  (triple gauge, all-outgoing convention)
#
# Ref: Peskin & Schroeder, Eq. (16.10) for triple gauge vertex
function _build_gauge_exchange(ch::TreeChannel, model::AbstractModel)
    (legL1, legL2), (legR1, legR2) = vertex_legs(ch)

    # Identify fermion pair vs boson pair
    if _both_fermion_legs(model, legL1, legL2)
        ferm1, ferm2 = legL1, legL2
        bos1, bos2 = legR1, legR2
    else
        ferm1, ferm2 = legR1, legR2
        bos1, bos2 = legL1, legL2
    end

    # Propagator Lorentz index (shared between fermion vertex and gauge vertex)
    rho = LorentzIndex(Symbol(:rho_, ch.channel), DimD())

    # Fermion chain: ū γ^ρ v
    chain = _fermion_line_chain(ferm1, ferm2, rho)

    # Triple gauge vertex in all-outgoing convention:
    # propagator carries q into the vertex → outgoing momentum is -q
    q = propagator_momentum(ch)
    neg_q = MomentumSum([(-c, m) for (c, m) in q.terms])
    mu1 = LorentzIndex(Symbol(:mu_, bos1.momentum.name), DimD())
    mu2 = LorentzIndex(Symbol(:mu_, bos2.momentum.name), DimD())
    vtx = triple_gauge_vertex(rho, mu1, mu2, neg_q, bos1.momentum, bos2.momentum)

    (chain, vtx)
end

# ---- Shared utilities ----

# Ref: Peskin & Schroeder, Section 4.6 (Feynman rules for fermions)
function _spinor_and_position(leg::ExternalLeg)
    p = leg.momentum
    m = leg.mass
    if leg.incoming && !leg.antiparticle
        (u(p, m), :right)     # incoming particle: u
    elseif leg.incoming && leg.antiparticle
        (vbar(p, m), :left)   # incoming antiparticle: vbar
    elseif !leg.incoming && !leg.antiparticle
        (ubar(p, m), :left)   # outgoing particle: ubar
    else
        (v(p, m), :right)     # outgoing antiparticle: v
    end
end
