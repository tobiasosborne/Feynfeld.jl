# Phase 18b-1c (bead feynfeld-gpi5): physical-antiparticle field labeling
# + spinor-dispatch fix.
#
# Step 1 (feynfeld-ghqq) wired ps1 through momentum routing. Step 2 fixes
# the field-label / spinor-dispatch convention so non-identity ps1
# emissions for Bhabha resolve to the correct physical spinors instead
# of crashing in walk_fermion_lines.
#
# Ground truth — qgraf "all incoming" convention:
#   refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:6961-6964 (qflow):
#     do i2=1,nleg
#       j3=flow(i1,invps1(i2))
#       if(i2.gt.inco)then; oflow(i1,i2)=-j3   ! outgoing — negated
#       else;               oflow(i1,i2)=j3    ! incoming — unchanged
#       end if
#     end do
#   refs/qgraf/ALGORITHM.md §1.2 (link array): conjugate fields are
#   represented as separate symbols (:e ↔ :e_bar).
#
#   In the all-incoming view, every leg flows INTO the diagram. The
#   per-slot field label flips between particle and antiparticle iff
#   the leg's physical direction is reversed: an outgoing physical
#   electron becomes a (would-be) incoming positron in this view.
#   Truth table for non-self-conjugate fields:
#     (anti=F, in=T)  in_e   → :e
#     (anti=T, in=T)  in_eb  → :e_bar
#     (anti=F, in=F)  out_e  → :e_bar    ← outgoing particle flips
#     (anti=T, in=F)  out_eb → :e        ← outgoing antiparticle flips
#   ⇒ slot label is :e_bar iff (antiparticle == incoming),
#     equivalently  iff (antiparticle XOR !incoming).
#
#   Spinor identity (u/v/ubar/vbar) must use the PHYSICAL flags, not the
#   qgraf-label-derived flags, because qgen's pmap reflects all-incoming
#   labels which disagree with physical antiparticle for outgoing legs.
#
# Cross-reference: refs/FeynCalc/.../ElAel-ElAel.m:77 uses
#   SetMandelstam[s, t, u, p1, p2, -k1, -k2, ...]
# i.e. physical outgoing momenta enter the Mandelstam invariants negated
# — the same all-incoming convention.

using Test
using Feynfeld
using Feynfeld: _qgraf_ext_labels
using Feynfeld.QgrafPort: _foreach_emission, emission_to_amplitude

@testset "Step 2: physical-antiparticle field labeling + dispatch fix" begin
    p1 = Momentum(:p1); p2 = Momentum(:p2); k1 = Momentum(:k1); k2 = Momentum(:k2)
    model = qed_model(m_e=:zero)

    @testset "_qgraf_ext_labels: XOR rule for all-incoming labels" begin
        in_e   = ExternalLeg(:e, p1, true,  false)   # → :e
        in_eb  = ExternalLeg(:e, p2, true,  true)    # → :e_bar
        out_e  = ExternalLeg(:e, k1, false, false)   # → :e_bar (flip)
        out_eb = ExternalLeg(:e, k2, false, true)    # → :e     (flip)
        @test _qgraf_ext_labels(model, [in_e, in_eb, out_e, out_eb]) ==
                [:e, :e_bar, :e_bar, :e]

        # ee→μμ (Phase 18a baseline) under XOR: matches the alternation
        # that diagram_gen's _expand_external_fields produces, so identity
        # ps1 is unaffected.
        in_mu  = ExternalLeg(:mu, k1, false, false)
        out_mub = ExternalLeg(:mu, k2, false, true)
        @test _qgraf_ext_labels(qed_model(m_e=:zero, m_mu=:zero),
                                [in_e, in_eb, in_mu, out_mub]) ==
                [:e, :e_bar, :mu_bar, :mu]

        # Self-conjugate (photon :gamma) is unchanged regardless of in/out.
        in_g  = ExternalLeg(:gamma, p1, true, false)
        out_g = ExternalLeg(:gamma, k1, false, false)
        @test _qgraf_ext_labels(model, [in_g, out_g]) == [:gamma, :gamma]
    end

    @testset "Bhabha: all 16 emissions resolve through emission_to_amplitude" begin
        # Pre-Step-2 (Session 31 diagnostic): of the 16 (ps1, pmap) tuples
        # qgen emits for Bhabha tree, 8 succeed (all routing s-channel by
        # coincidence of identity-ps1) and 8 throw "multiple plain-spinor
        # legs" / "multiple bar-spinor legs" from walk_fermion_lines —
        # because _spinor_dispatch reads the qgraf-label-derived anti
        # flag, which mis-classifies an outgoing electron at a slot whose
        # qgraf label is :e_bar. After Step 2 both the field labels (XOR
        # rule via _qgraf_ext_labels) AND the dispatch (phys_anti threaded
        # through emission_to_amplitude / build_externals) are physical-
        # convention-aware, so all 16 emissions yield AmplitudeBundles.
        #
        # solve_tree_pipeline currently masks the failure via is_emission_
        # canonical (the Strategy-C filter, bead vjw9), which keeps only
        # 1 of the 2 orbits — so this test bypasses the filter and walks
        # _foreach_emission directly.
        in_e   = ExternalLeg(:e, p1, true,  false)
        in_eb  = ExternalLeg(:e, p2, true,  true)
        out_e  = ExternalLeg(:e, k1, false, false)
        out_eb = ExternalLeg(:e, k2, false, true)
        legs = [in_e, in_eb, out_e, out_eb]
        physical_moms = [leg.momentum for leg in legs]
        phys_anti = Bool[leg.antiparticle for leg in legs]
        ext_exp   = _qgraf_ext_labels(model, legs)
        in_raw    = [in_e.field_name, in_eb.field_name]
        out_raw   = [out_e.field_name, out_eb.field_name]

        success      = Ref(0)
        denom_strs   = Set{String}()
        _foreach_emission(model, in_raw, out_raw; loops=0,
                            ext_exp=ext_exp) do state, labels, ps1, pmap
            bundle = emission_to_amplitude(state, labels, ps1, pmap, model;
                                             physical_moms=physical_moms,
                                             n_inco=2, phys_anti=phys_anti)
            success[] += 1
            for d in bundle.denoms
                push!(denom_strs, string(d))
            end
        end
        @test success[] == 16

        # Acceptance #4 (bead gpi5): post-Step-2, the 16 emissions surface
        # both s-channel (p1+p2)² and t-channel (p1-k1)²/(p2-k2)² propagator
        # denominators. Pre-Step-2: every emission routed s only (8 s + 8
        # crashes). Pre-Step-1: ALL emissions routed (p1+p2)² because ps1
        # was ignored. So distinct-denom count > 1 is the load-bearing
        # assertion that ps1-threaded routing + physical-anti dispatch
        # together expose both channels.
        @test length(denom_strs) > 1
        # Spot-check t-channel presence: at least one denom involves
        # (p1, k1) (or its conjugate-pair partner (p2, k2)) — evidence
        # that the t-channel topology v5=(p1,k1), v6=(p2,k2) emits.
        has_t = any(s -> (occursin("p1", s) && occursin("k1", s)) ||
                          (occursin("p2", s) && occursin("k2", s)),
                     denom_strs)
        @test has_t
    end

    @testset "Phase 18a regression: ee→μμ pipeline ≡ handbuilt symbolically" begin
        # Identity-ps1 path unchanged by Step 2 — XOR rule and alternation
        # agree for ee→μμ since each species occurs once in + once out.
        in_e   = ExternalLeg(:e,  p1, true,  false)
        in_eb  = ExternalLeg(:e,  p2, true,  true)
        out_mu = ExternalLeg(:mu, k1, false, false)
        out_mub= ExternalLeg(:mu, k2, false, true)
        prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                                    [in_e, in_eb], [out_mu, out_mub], 10.0)
        r_pipeline  = solve_tree_pipeline(prob)
        r_handbuilt = solve_tree(prob)
        @test r_pipeline.amplitude_squared == r_handbuilt.amplitude_squared
    end
end
