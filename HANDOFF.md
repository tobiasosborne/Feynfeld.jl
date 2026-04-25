# HANDOFF — 2026-04-25 (Session 32: vjw9 Step 2 landed)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE (Session 32 updates)

1. Read `CLAUDE.md` first. Then this Session 32 block. Then Session 31 for the
   4-step Bhabha dependency chain.
2. **Step 2 (`feynfeld-gpi5`) landed and closed.** Bhabha emissions now route
   correctly: 16 (ps1, pmap) tuples all resolve through `emission_to_amplitude`
   without crashing, and the denom set partitions into both s `(p1+p2)²` and
   t `(p1-k1)²` channels — acceptance criterion #4 met.
3. **Default next work: `feynfeld-vjw9` (Step 4, orbit dedup).** With Step 2
   in place, the 16 emissions split into a 2-orbit structure that is now
   visible to the rest of the pipeline. The original Strategy A/B/C options
   apply afresh; Session 31's open question on re-running `same_emission_orbit`
   under Step-2 routing is now actionable. `solve_tree_pipeline(Bhabha)` still
   returns `n_emissions=1` because `is_emission_canonical` drops the t-orbit's
   canonical rep — that's the remaining vjw9 work.
4. **Phase 18b-3 (`feynfeld-a7f2`, Step 3)** stays NOT on Bhabha's critical
   path. Independent; required for Compton's internal fermion propagator.

## SESSION 32 TIMELINE

1. Read HANDOFF.md + Feynfeld_PRD.md. Summarised project rules + Step 2 plan.
   Tobias: "go for it. (As long as it conforms to eg qgraf ground truth)".
2. Dispatched 4 parallel Sonnet Explore agents for ground-truth verification:
   Step 1 wiring as it stands, spinor-dispatch + field-expansion plumbing,
   Bhabha test ground truth (FeynCalc + handbuilt), qgraf ps1/pmap convention
   incl. f08:6961-6964 cross-check. Synthesised into a 4-file Step 2 plan.
3. Tobias approved. Claimed `feynfeld-gpi5`. Baselined RED test
   `test/v2/qgraf/test_phase18b1_multi_orbit.jl` (2 pass + 2 broken, 18s).
4. Wrote new RED test file `test/v2/qgraf/test_phase18b1_step2_dispatch.jl`
   with 3 testsets: XOR-rule helper unit, all-16-Bhabha-emissions resolve
   through direct `_foreach_emission` walk, Phase 18a regression. Verified
   RED via 4 errors (`_qgraf_ext_labels` undefined). One assertion (test 2)
   initially passed accidentally because `is_emission_canonical` filters
   `solve_tree_pipeline` down to 1 emission — replaced with direct enum.
5. Implemented Step 2 across 5 files (~115 LOC diff):
   - `src/v2/cross_section.jl` — added `_qgraf_ext_labels(model, legs)` helper
     with the XOR rule + ground-truth citations; `solve_tree_pipeline` now
     builds `phys_anti` and `ext_exp`, threads both through.
   - `src/v2/qgraf/audition.jl::_foreach_emission` — added
     `ext_exp::Union{Nothing,Vector{Symbol}}=nothing` kwarg; bypasses
     `_expand_external_fields` when supplied.
   - `src/v2/qgraf/emission_amplitude.jl::emission_to_amplitude` — added
     `phys_anti::Union{Nothing,Vector{Bool}}=nothing` kwarg, threaded to
     `build_externals` + `walk_fermion_lines`.
   - `src/v2/qgraf/vertex_assemble.jl::build_externals` — added `phys_anti`
     kwarg; in the loop, `anti = phys_anti===nothing ? _is_antiparticle_field(field)
     : phys_anti[phys_idx]`. Default preserves pre-Step-2 behaviour for
     test/older callers.
   - `src/v2/qgraf/fermion_line.jl::walk_fermion_lines` — added `phys_anti`
     kwarg, forwarded.
6. Re-ran the new test: 5/5 green at first cut. Strengthened with two more
   assertions (denom-set size > 1 and t-channel-marker check) after a quick
   `/tmp/step2_denom_check.jl` diagnostic showed 16 emissions and 4 distinct
   denom string forms covering both s `(p1+p2)²` and t `(p1-k1)²`/`(p2-k2)²`.
   Final test count: 7 assertions, all green.
7. Multi-orbit RED test still 2 pass + 2 broken — correctly, those are vjw9's
   gates which Step 4 unblocks.
8. Full suite: 1332 pass + 6 broken (was 1325+6 baseline). +7 new, 0 regressions.
9. Reviewer agent (general-purpose, given full diff context + ground-truth
   citations to verify): NO BLOCKERS. 3 deferred concerns:
   - C1 (semantic overload) — `ExternalFactor.antiparticle` field semantics
     shifted from label-derived to physical-when-supplied. Worth a follow-up
     to add an explicit `phys_antiparticle::Bool` field or rename. Not
     blocking; the docstring acknowledges the duality.
   - C2 (idiom) — `Union{Nothing, Vector{Bool}}` 3 layers deep is awkward.
     Could collapse if pipeline always synthesises `phys_anti`. Style.
   - C3 (drift risk) — `count_dedup_burnside_qg21` and `count_dedup_prefilter`
     in `audition.jl` still use `_expand_external_fields`. Counts depend only
     on the multiset, so this is correct, but next reader could propagate
     the Bhabha mistake back. **TODO comment added** at audition.jl:316.
10. Picked up nit (error message wording at audition.jl:212).
11. Closed `feynfeld-gpi5`, updated `feynfeld-vjw9` notes with Step 2 status
    + Step 4 re-evaluation suggestion.
12. This HANDOFF. `bd export -o .beads/issues.jsonl`. Commit + push.

## SESSION 32 ACCOMPLISHMENTS

- **Step 2 of the 4-step Bhabha unblock chain landed cleanly.** No regressions,
  full reviewer pass, ground-truth citations on every new docstring.
- **Bhabha emissions now physically routable end-to-end up to orbit dedup.**
  s + t channels both surface in the emission denom signatures. The remaining
  Bhabha blocker is genuinely just orbit selection (vjw9), not anything
  upstream in the pipeline.
- **Reviewer-flagged follow-ups documented but deferred.** None gate vjw9.

## SESSION 32 OPEN QUESTIONS (for next agent)

- **Step 4 strategy choice.** With Step 2 in place, `same_emission_orbit`
  (audition.jl, Session 29) on Bhabha may now give the correct 2-orbit
  partition (it gave 4 under the broken routing). Re-run before re-attempting
  Strategy A vs C. If 2 orbits, Strategy C (orbit-grouping with one rep per
  orbit) is the simplest path; if not, Strategy B (Burnside-all + relabeling
  in spin_sum_interference, bead `feynfeld-rj1l`) remains the alternative.
- **C1 follow-up bead?** Worth filing a small cleanup bead for after the
  Bhabha chain closes: add explicit `ExternalFactor.phys_antiparticle::Bool`
  + remove the label-derived fallback in `build_externals`. Defer until
  Steps 4/5 land.

---

# HANDOFF — 2026-04-21 (Session 31: vjw9 scope correction + Step 1 landed)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE (Session 31 updates)

1. Read `CLAUDE.md` first. Then this Session 31 block. Then Session 30 for vjw9 context.
2. **Session 30's diagnosis was incomplete.** The ps1-threading fix it identified is one
   of FOUR interlocking steps to unblock Bhabha. Step 1 landed this session; Steps 2-5
   remain. The test file `test/v2/qgraf/test_phase18b1_multi_orbit.jl` still shows
   2 `@test_broken` — that's expected; vjw9 isn't closed.
3. **The real 4-step dependency chain** for Bhabha (each a discrete, testable unit):
   - **Step 1 — ps1-threading + ext_signs outgoing-negation** [`feynfeld-ghqq`, DONE].
     Thread ps1 through emission_to_amplitude → build_externals → walk_fermion_lines;
     add ext_signs kwarg to route_momenta for the qgraf all-incoming convention.
     Matches qgraf Fortran `qflow` f08:6961-6964 (verified). Identity default preserves
     existing callers. Single-bundle Bhabha pipeline now routes s-channel correctly
     (produces `T_ss` exactly). ee→μμ regression green. 1325 tests pass, 0 regressions.
   - **Step 2 — Physical-antiparticle field labeling + dispatch fix** [`feynfeld-gpi5`, open,
     blocked by Step 1's commit but now unblocked]. `cross_section.jl:163` drops
     `leg.antiparticle`; `_expand_external_fields` alternation happens to work for
     identity-ps1 by coincidence but corrupts t-channel orbit. `_spinor_dispatch` uses
     field-label-derived anti flag, giving `v(k1)` for outgoing electron instead of
     `ubar(k1)`. Fix: XOR rule for ext_raw labels + thread phys_antiparticle through
     dispatch. ~30-50 LOC. Critical insight: **with Step 2 done, Bhabha's t-channel
     fermion lines stay 1-vertex** — Phase 18b-3 (`feynfeld-a7f2`) is NOT on Bhabha's
     critical path (it's for Compton's internal fermion propagator).
   - **Step 3 — walk_fermion_lines multi-vertex** [`feynfeld-a7f2`, Phase 18b-3].
     Required for Compton, NOT Bhabha. Independent of vjw9.
   - **Step 4 — Orbit dedup** [`feynfeld-vjw9`, scope corrected]. With Step 2 done,
     16 emissions will partition into 2 orbits (8 s + 8 t, not 8 s + 8 u as under
     current convention). Original Options A/B/C apply. May need fresh empirical
     test post-Step 2 to pick.
   - **Step 5 — Interference label matching** [`feynfeld-rj1l`]. Scope depends on
     Step 4. May be trivial if orbit reps have consistent bar_mom labels.
4. **Default next work**: Step 2 (`feynfeld-gpi5`). Self-contained, ~30-50 LOC, has
   clear empirical test. Unblocks vjw9 if combined with a Step 4 decision.

## SESSION 31 TIMELINE

1. Read HANDOFF.md + Feynfeld_PRD.md. Summarized project rules for Tobias.
2. Dispatched 5 parallel Sonnet Explore agents to map the codebase for vjw9:
   call-chain audit, qgen/ps1/pmap semantics, Bhabha ground truth, interference plumbing,
   Julia idioms + existing callers. Synthesized into concrete fix plan.
3. Dispatched 4 more parallel agents for ground-truth verification per Rule 3:
   qgraf Fortran momentum convention, FeynCalc Bhabha handling, route_momenta audit,
   QEDcore.jl ecosystem sweep. Verified the outgoing-negation claim against qgraf's
   own `qflow` sign-flip at f08:6961-6964 (not just the docstring).
4. Tobias approved the plan: \"ok go for it. red green tdd\".
5. Baselined RED test: 2 pass + 2 @test_broken (13.2s). Implemented Step 1:
   - `momentum.jl::route_momenta` — `ext_signs::Vector{Int}=fill(1, n_ext)` kwarg;
     applied as coefficient multiplier in terms loop (line 220).
   - `vertex_assemble.jl::build_externals` — `ps1::AbstractVector{<:Integer}=1:n_ext`
     kwarg; uses `phys_idx=Int(ps1[i])` for both `physical_moms[phys_idx]` and
     `incoming = phys_idx <= n_inco`. leg_idx remains slot index.
   - `fermion_line.jl::walk_fermion_lines` — same ps1 kwarg; threads to internal
     `build_externals` call.
   - `emission_amplitude.jl::emission_to_amplitude` — builds `qgraf_ext_moms[i] =
     physical_moms[Int(ps1[i])]` + `ext_signs[i] = phys_idx > n_inco ? -1 : 1`;
     passes `ps1=ps1` kwarg to build_externals and walk_fermion_lines. Validates
     ps1 length.
6. Ran RED test → STILL 2 pass + 2 @test_broken. Diagnosed: pipeline still reports
   `n_emissions=1` because `is_emission_canonical` filters. Ran diagnostic: single
   bundle produces **exactly T_ss** (verified via `evaluate_sp` against handbuilt
   trace). Momentum routing works — filter is the remaining blocker.
7. Ran full suite (5min): 1325 pass + 6 broken, no regressions.
8. Enumerated all 16 Bhabha emissions: 8 succeed (all s-channel), 8 crash in
   walk_fermion_lines with \"multiple plain-spinor legs\" / \"multiple bar-spinor legs\".
9. Dug into crash root cause: certain ps1 values produce vertices with both slots
   at the same chain position (both :right or both :left). Initially thought this
   meant Phase 18b-3 multi-vertex fermion lines (feynfeld-a7f2) was on critical path.
10. Diagnostic test: forced physical ext_exp=[:e,:e_bar,:e_bar,:e] (bypassing
    _expand_external_fields alternation). Crashes still happened because
    `_spinor_dispatch` uses field-label-derived anti. Traced that to `cross_section.jl:163`
    dropping is_antiparticle from in_fields — the real root cause is Step 2, not Step 3.
11. Corrected the plan with Tobias's approval: 4-step dependency chain with Step 2
    (labeling + dispatch) between Step 1 and Step 4.
12. Reviewer agent on Step 1 changes: no bugs found; 3 acceptable notes (ext_signs
    allocation in hot path, AbstractVector{<:Integer} polymorphism, no non-identity
    unit test). Committed.
13. Filed new beads: `feynfeld-ghqq` (Step 1, close now) and `feynfeld-gpi5` (Step 2).
    Wired dependencies: gpi5 depends on ghqq, vjw9 depends on gpi5. Updated vjw9
    bead notes with Session 31 findings.
14. Updated HANDOFF.md (this block).
15. `bd dolt push` + `git push`.

## SESSION 31 ACCOMPLISHMENTS

- **Step 1 landed, foundational and correct.** 4 files touched (~40 LOC changed).
  Cross-validated against qgraf Fortran source (f08:6961-6964), FeynCalc's
  `SetMandelstam[..., -k1, -k2]` convention, and `qdis_fermion_sign`'s own ps1 usage.
- **Scope correction for vjw9.** The 4-step dependency chain is now documented in
  the bead, in beads for Steps 1/2, and in this HANDOFF. Future agents won't
  chase the wrong fix.
- **Bhabha critical path clarified.** Phase 18b-3 (multi-vertex fermion lines) is
  NOT required for Bhabha; only for Compton. Significant scope reduction for Bhabha.
- **Empirical verification pattern established.** /tmp/diag_* scripts (not committed)
  let me test hypotheses without modifying source. Pattern: directly invoke
  `_foreach_emission` with manual ext_exp override, catch errors, dump denom
  signatures. Reproducible by next session if needed.

## SESSION 31 OPEN QUESTIONS (for next agent)

- **Step 2's exact shape.** The bead gpi5 outlines the fix but leaves some details
  unspecified: should `_expand_external_fields` be fixed to not alternate when input
  is already correct, or bypassed entirely? Does ExternalFactor need both `antiparticle`
  (phys) and a derived `is_bar_field` (from pmap label), or just phys_antiparticle?
  Needs a focused design session with the reviewer at the end.
- **Step 4 re-evaluation.** The Strategy-C bug description in vjw9 (canonical rep
  may be qgen-invalid) was under the OLD convention. After Step 2, the orbit
  structure may look different — possibly `is_emission_canonical` works fine without
  any extra fix. Re-run `count_dedup_canonical` after Step 2 before picking strategy.
- **HANDOFF's \"multi-vertex walk_fermion_lines\" scope.** Currently a7f2's description
  says \"walk the fermion line across internal fermion propagators\" — correct for
  Compton but its dependency on h3pb (fermion propagator composite momentum) is
  right. Bhabha is NOT a use case for a7f2; scope note added to vjw9 but a7f2's
  description is fine as-is.

## START HERE (Session 30 updates)

1. Read `CLAUDE.md` first. Then this Session 30 block. Then Session 29.
2. **`feynfeld-vjw9` blocker is NOT an orbit-dedup bug.** The original A/B/C options
   in the bead (and Session 29's attempts) targeted `audition.jl`. Session 30
   diagnostic shows the real bug is upstream in `src/v2/qgraf/emission_amplitude.jl:64-68`:
   `qgraf_ext_moms[i] = physical_moms[i]` ignores `ps1`, so every Bhabha emission
   routes s-channel `(p1+p2)²` regardless of which ps1 permutation it came from.
   Pipeline produces `|M|² == T_ss` exactly — zero t-channel contribution.
3. Test suite unchanged: 2 pass + 2 `@test_broken` in
   `test/v2/qgraf/test_phase18b1_multi_orbit.jl`. Phase 18a regression green.
   **No source code modified this session.**
4. Next agent should read the Session 30 block below in full, then either
   (a) do the architectural ps1-threading fix (non-trivial, multi-file), or
   (b) pick a different ready bead and defer vjw9 until we've discussed the
   architecture with Tobias.

## SESSION 30 TIMELINE — diagnostic-only session, no code changes

1. Tobias: "read HANDOFF.md, Feynfeld_PRD.md … then get to work with whatever
   you think best". Dispatched three parallel read-only Sonnet Explore agents
   on the vjw9 domain (audition.jl logic, Bhabha ground truth, burnside_combine
   + interference plumbing).
2. Initial hypothesis: Option A (lex-next-smallest fallback) — ~15 LOC. Agent-1
   recommended orbit-signature dedup. But the bead notes said Session 29 tried
   exactly that and it failed empirically ("theoretical cocycle identity doesn't
   hold through sorting"). Correct — I needed to see real data.
3. Baselined the RED tests (2 pass + 2 `@test_broken`, 13s).
4. Wrote six diagnostic scripts into `/tmp/diag_*.jl` (not committed):
   - `diag_bhabha_orbits.jl`: all 16 emissions + same_emission_orbit partition.
   - `diag_orbit_via_diagsig.jl`: `_diagram_sig`-based partition (gave 5 classes).
   - `diag_full_orbit.jl`: full G-orbit of each emission, showing X is not G-closed.
   - `diag_denoms.jl`: denom signature per canonical emission.
   - `diag_permute_test.jl`: route_momenta with ps1-permuted moms (proved t-channel
     EMERGES once you permute physical_moms by ps1).
   - `diag_all16_denoms.jl`: all 16 denoms under ps1-permutation (4 classes: s,
     k1+k2, t, p2+k2 — 4 not 2 because route_momenta doesn't negate outgoing).
   - `diag_correct_orbits.jl`: redo with the CORRECT `[:e,:e]` convention (see §4
     below) — gave 5 orbit-sig groups all with denom (p1+p2)².
5. Followed a rabbit hole for most of the session using
   `in_fields=[:e,:e_bar],[:e,:e_bar]` — the field-name expansion path that
   DOESN'T match what `solve_tree_pipeline` actually does. Caught the error only
   at the end when investigating why my in-callback diagnostic said 3 canonical
   but the actual pipeline reports `n_emissions=1`.
6. Updated bead `feynfeld-vjw9` with corrected findings (twice — first note
   was based on wrong convention, second note supersedes it).
7. Did NOT change any source files. Committed only `.beads/interactions.jsonl`
   with the Session 30 narrative. Pushed.

## SESSION 30 ACCOMPLISHMENTS — diagnostic understanding, zero code

- **Located the real bug.** `emission_amplitude.jl:64-68` does not thread `ps1`
  through `physical_moms`. The `for i in 1:n_ext; qgraf_ext_moms[i] = physical_moms[i]; end`
  loop is a no-op rename. Under this, every ps1 permutation routes the same
  slot-indexed physical momenta → same propagator momentum → all s-channel for
  Bhabha.  Proved empirically: with a manual `permuted[i] = physical_moms[ps1[i]]`
  in a standalone script, t-channel denom `(p1+k1)` appears for the non-identity
  ps1 permutations.
- **Documented why orbit-dedup in audition.jl cannot fix vjw9.** X=16 emissions
  is not G-closed. `same_emission_orbit` splits 16 into {8,4,2,2}; `_diagram_sig`
  lex-min gives 5 classes; both conventions' orbit-sig groupings give 5. Burnside
  arithmetic `16·|Stab|/|G| = 16·1/8 = 2` equals the physical channel count by
  design, not by partitioning X. Any audition.jl-level "pick 2 canonical reps"
  is a dead end because the 16 bundles currently route identical s-channel
  momenta — you cannot distinguish s from t by inspecting bundles alone.
- **Corrected stale bead description.** Prior sessions (22, 27, 29) framed vjw9
  as a canonical-rep-may-be-qgen-invalid Strategy-C bug. The Strategy-C
  description is still accurate (canonical filter gives n_emissions=1, not 2),
  but the fix is NOT in `is_emission_canonical`; it is upstream in
  `emission_to_amplitude`.

## SESSION 30 FINDINGS — the architectural bug in detail

### 1. Topology and emission count

qg21 generates **one** topology state for Bhabha tree (not two). Adjacency:
v1, v2 → v5; v3, v4 → v6; single v5-v6 photon edge. All 16 emissions share this
vmap. `_foreach_emission` yields 16 `(ps1, pmap)` tuples per this state.
`|autos| = 8`, every emission has `|Stab| = 1`, so Burnside arithmetic gives
`16/8 = 2`.

### 2. The convention gotcha that cost most of the session

`solve_tree_pipeline` (cross_section.jl:163-164) does:
```julia
in_fields  = [leg.field_name for leg in prob.incoming]
out_fields = [leg.field_name for leg in prob.outgoing]
```

`ExternalLeg.field_name` is `:e` regardless of `is_antiparticle`. So for a
`CrossSectionProblem` set up as standard Bhabha (`[in_e, in_eb], [out_e, out_eb]`),
`solve_tree_pipeline` calls `_foreach_emission` with `in=[:e,:e], out=[:e,:e]`,
NOT `[:e,:e_bar],[:e,:e_bar]`. `_expand_external_fields` then alternates
conjugates: `[:e,:e,:e,:e]` → `ext_exp=[:e,:e_bar,:e,:e_bar]`.  This is the
"all particles, antiparticle-ness inferred" qgraf convention. Valid, but
different from `ext_raw=[:e,:e_bar,:e,:e_bar]` which expands to
`ext_exp=[:e,:e_bar,:e_bar,:e]`.

Both conventions produce 16 emissions and Burnside count 2 — but the SPECIFIC
ps1/pmap tuples differ. When I analysed `[:e,:e_bar]` manually for most of
the session, the canonical count was 3, misleading me into thinking
`is_emission_canonical` was over-accepting. Under the pipeline's actual
convention `[:e,:e]`, canonical count is 1 (the original Strategy-C under-count).

**Convention used in production code and reference: `[:e,:e],[:e,:e]`.**

### 3. Why all bundles currently produce s-channel

`emission_to_amplitude` (emission_amplitude.jl:54-105) receives
`physical_moms = [p1, p2, k1, k2]` in slot order and `ps1` from qgen. Lines
64-68:
```julia
qgraf_ext_moms = Vector{Momentum}(undef, n_ext)
for i in 1:n_ext
    qgraf_ext_moms[i] = physical_moms[i]  # no ps1 permutation
end
```
The comment says "all incoming by negating outgoing" but the code neither
permutes by ps1 nor negates outgoing. `route_momenta` then sees
`[p1,p2,k1,k2]` in slot order regardless of ps1, routes propagator
`v5→v6` = `p1+p2` always.

Diagnostic confirmation (`/tmp/diag_permute_test.jl`): when I manually built
`permuted[i] = physical_moms[ps1[i]]` and called `route_momenta` on it,
non-identity ps1 permutations gave `(p1+k1)²` denoms — the t-channel.

### 4. Why orbit-based dedup can't rescue this

The current code gives ALL 16 emissions the same denom (p1+p2)². Dedup at the
bundle level by denom signature yields 1 class, not 2. Dedup by orbit
equivalence at the `(ps1, pmap)` level yields 5 classes (visible sizes
{6,4,2,2,2} across the `[:e,:e]` convention). Neither matches the physical 2.

X=16 is not G-closed: for emission e, `|Orbit_G(e)| = |G|/|Stab(e)| = 8` but
many of those 8 images are not themselves valid qgen emissions (qgen picks
one pmap per `(state, ps1)`; group-action images of a chosen pmap are
mathematically valid but qgen never emits them).

Even `count_dedup_canonical` (Strategy B) gives 1, confirming that picking
"the canonical orbit rep per lex-min" doesn't yield 2 reps — one orbit's
lex-min lies outside X.

Any fix that stays inside `audition.jl` is chasing ghosts. The t-channel
must be SURFACED at the `emission_to_amplitude` level by threading `ps1`
through momentum routing.

### 5. The actual fix plan (NOT implemented Session 30)

Step 1. Patch `src/v2/qgraf/emission_amplitude.jl:54` `emission_to_amplitude`:
```julia
qgraf_ext_moms = Vector{Momentum}(undef, n_ext)
for i in 1:n_ext
    qgraf_ext_moms[i] = physical_moms[Int(ps1[i])]
end
```
And pass `ps1` through to `build_externals` and `walk_fermion_lines`.

Step 2. Add `ps1` kwarg (default `1:n_ext` — identity) to
`src/v2/qgraf/vertex_assemble.jl::build_externals`. Inside:
```julia
phys_idx = Int(ps1[i])
field    = pmap[i, 1]
mom      = physical_moms[phys_idx]
incoming = phys_idx <= n_inco  # was: i <= n_inco
```
This moves incoming-flag from slot-based to ps1-remapped. Existing callers
(two test files: `test_fermion_line.jl`, `test_external_assemble.jl`) don't
pass `ps1`, so they default to identity and behave unchanged.

Step 3. Mirror ps1 kwarg through `src/v2/qgraf/fermion_line.jl::walk_fermion_lines`
(line 45) to `build_externals`.

Step 4. Run `test/v2/runtests.jl` — identity ps1 for ee→μμ means Phase 18a
regression should be no-op. The 30 qgraf tests should be unaffected.

Step 5. Bhabha will then emit bundles with TWO distinct denom classes
(s: p1+p2 and t: p1+k1 up to conjugate-pair shifts). Under the original
`[:e,:e_bar]`-expansion convention my diagnostic showed 4 denom classes
(s, -s, t, -t). Under the correct `[:e,:e]` convention this needs re-running
post-fix. Dedup to 2 reps via either (a) denom-signature grouping with
sign-normalised momenta, or (b) diagonal-trace AlgSum equality. Then
`burnside_combine` with weights=1.

Step 6. Interference between s-rep and t-rep crosses the
`_find_line_by_bar_mom` convention boundary that Session 29 hit. If Step 5
produces reps whose fermion-line momenta line up, it may Just Work; if not,
fold in the Option-B momentum-relabeling fix (bead `feynfeld-rj1l`).

Step 7. Validate symbolically against the handbuilt
`T_tt + T_ss − 2·T_int` reference in `test/v2/test_bhabha.jl`. Unlocks the
two `@test_broken` assertions and closes vjw9.

### 6. Open risks in the fix plan

- **Spinor convention under ps1 permutation.** If ps1 places an
  originally-outgoing leg (`ps1[i] > n_inco`) at an incoming slot, the
  field at that slot (from `pmap[i,1]`) may have the "wrong" particle/
  antiparticle sign relative to physical Bhabha. With conjugate-expanded
  `ext_exp`, the field is already flipped for alternating same-species
  slots. Net spinor construction will need to be verified emission by
  emission against handbuilt reference; this is where the change can
  silently break |M|² without throwing.
- **Fermion sign.** `_emission_fermion_sign` uses ps1 internally;
  unclear whether its current convention is consistent with the proposed
  ps1-threaded physical_moms. Likely needs cross-checking.
- **No outgoing negation.** The pre-existing "negate outgoing" comment at
  line 47-48 is still a no-op after the fix. `route_momenta` may rely on
  all-incoming convention internally; if so, post-ps1-permutation slots
  need sign flip for legs where `ps1[i] > n_inco`. Or route_momenta may
  not care (|M|² is even in momenta). Untested.

### 7. Diagnostic scripts — kept in /tmp, NOT committed

| Script | Purpose |
|--------|---------|
| `/tmp/diag_bhabha_orbits.jl` | All 16 emissions with ps1/pmap/canonical/stab dump |
| `/tmp/diag_orbit_via_diagsig.jl` | `_diagram_sig`-based orbit partition (5 classes) |
| `/tmp/diag_full_orbit.jl` | Full G-orbit under auto action; shows X not G-closed |
| `/tmp/diag_denoms.jl` | Denom signature per canonical emission (all s-channel) |
| `/tmp/diag_permute_test.jl` | route_momenta with ps1-permuted moms — proves t-channel emerges |
| `/tmp/diag_all16_denoms.jl` | All 16 denoms under ps1-permutation (4 classes) |
| `/tmp/diag_current_amp.jl` | Pipeline |M|² == T_ss (1 emission only, s-channel) |
| `/tmp/diag_correct_orbits.jl` | Correct `[:e,:e]` convention orbit analysis |
| `/tmp/diag_state_count.jl` | Confirms 1 qg21 state for Bhabha |
| `/tmp/diag_inside_callback.jl` | Correct 3-canonical count using `[:e,:e_bar]` (wrong convention; kept for comparison) |
| `/tmp/diag_pipeline_trace.jl` | Bundle-accumulation trace under try/catch |

These exist only in /tmp (ephemeral). If the next session wants to reproduce,
they'll need to recreate. Most are 20-40 LOC; regenerating is cheap.

## SESSION 30 HANDOFF — what the next agent should do

1. **Read the Session 30 findings block above in full.** The bug is
   architectural, not in audition.jl. Understand that before touching code.
2. **Decide the path with Tobias.** Options:
   - **A (fix vjw9 now).** Work through the 7-step plan in §5. ~200-300 LOC
     across 3 files (emission_amplitude.jl, vertex_assemble.jl,
     fermion_line.jl) plus test verification. Risk: spinor convention
     under permutation — verify step by step against handbuilt.
   - **B (defer vjw9, work elsewhere).** Pick from the ready queue:
     Move 1.4f @inferred regression (`feynfeld-4206`), 1.6 tutorial (`lj1`),
     MUnit batches (`q6m`, `s4p`, `8qe`, `n01`, `iaz`), Eps contraction
     (`qu1`). Most Phase 18b beads (`h3pb`, `m4o8`, etc.) are gated on
     vjw9; they remain blocked.
   - **C (discuss architecture with Tobias first).** The ps1-threading
     gap has been latent since Phase 18a landed. If Tobias has a preferred
     approach (e.g. make n_inco a per-slot flag; rework qgen's pmap
     emission to be channel-explicit; something else), the fix should
     align with that direction.
3. **Session close protocol.** As always: `bd dolt push` + `git push`.

---

1. Read `CLAUDE.md` first — the v1 paragraph was rewritten Session 29 to reflect
   that `src/algebra/`, `src/integrals/`, `src/Feynfeld.jl`, `test/algebra/`,
   `test/integrals/`, `test/runtests.jl`, `test/test_ee_mumu.jl`, and the empty
   `src/{model,rules,diagrams,evaluate}/` + `test/{model,rules,diagrams,evaluate,references}/`
   scaffolds are **deleted**, not "will be deleted".
2. `test/v2/runtests.jl` is the canonical aggregate entry point — single-process,
   605 assertions, ~5 min. Use this, NOT `grind/run_v2_tests.sh` (fork-per-file).
3. Phase 18b blocker unchanged: `feynfeld-vjw9` (orbit-rep dedup, `audition.jl:69`).
   Pick this up next unless Tobias redirects.

## SESSION 28.5 CRASH SUMMARY (what happened, for forensics)

An unnamed Session 28.5 agent was executing bead `feynfeld-tzgc` (Move 1.1,
v1 deletion per Session 28 stocktake §10). Sequence before crash:
1. Created beads `tzgc` / `8cgv` / `4206`; claimed `tzgc`.
2. Pre-flight grep: zero cross-deps from `src/v2/` or `test/v2/` on `src/algebra/`
   or `src/integrals/`. Confirmed safe.
3. Baseline snapshot: ran `bash grind/run_v2_tests.sh > /tmp/pre_delete.log`.
4. `git rm` executed: the full v1 file set staged for deletion.
5. Post-delete regression: re-ran `bash grind/run_v2_tests.sh` — **exit 137
   (OOM SIGKILL) from the kernel**. WSL2 crashed.

**Root cause assessment** (Session 29, Tobias): concurrent load from a separate
project on the same WSL2 instance, not `grind/run_v2_tests.sh` per se. But the
script DOES fork N Julia processes (~300-500 MB each) and compounds pressure.
**Recommendation:** prefer `test/v2/runtests.jl` for aggregate runs — it loads
`FeynfeldX` once and reuses the JIT across all files, bounded memory.

Staged deletions survived the crash (git index is durable). `/tmp/pre_delete.log`
did not.

## SESSION 29 TIMELINE — pick up the staged deletion, verify, commit

1. Forensic reconstruction: dispatched a Sonnet subagent to read the crashed
   session JSONL (`~/.claude/projects/.../f136a661-....jsonl`, ~964 KB) and
   summarise the trajectory + OOM root cause + what was staged and why.
   Kept the raw JSONL out of the main context.
2. Confirmed `feynfeld-tzgc` still in-progress and matched the stocktake plan.
   Confirmed staged deletions correspond exactly to stocktake §10.
3. Verified acceptance grep: zero references to `src/algebra` / `src/integrals`
   in `src/v2/` or `test/v2/`. Only docs (CLAUDE.md, HANDOFF.md, Feynfeld_PRD.md)
   mention them historically — correct.
4. Ran `julia --project=. test/v2/runtests.jl` in a single process.
   **Result: 605 / 605 pass, 5m04.6s.** All green post-deletion.
5. Updated CLAUDE.md §"Active code" — v1 paragraph rewritten to past tense
   with Session 29 date and the exact deleted paths. Updated test/LOC numbers
   (301 → 605, 28 files → 69 files, 3,400 LOC → 10,300 LOC) to match the
   stocktake. Kept the "do not resurrect v1 patterns" directive.
6. Wrote this Session 29 block.
7. Committed staged deletions + doc updates as a single commit for bisectability.
8. `bd close feynfeld-tzgc` (acceptance met).
9. `bd dolt push` + `git push`.

## SESSION 29 ACCOMPLISHMENTS

- v1 deletion landed cleanly: ~5,400 LOC removed, 605/605 v2 tests green (commit `db3399d`).
- `feynfeld-tzgc` closed. Staged work from the crashed Session 28.5 recovered with zero loss.
- CLAUDE.md §"Active code" resynced with stocktake reality.
- Forensic reconstruction process documented (Sonnet subagent on crashed JSONL) so future
  agents have a pattern for WSL-crash recovery.

### Session 29 addendum (later in same session): Moves 1.2 + 1.3 landed

- **Move 1.2 (`feynfeld-qyu`, FeynfeldX → Feynfeld rename).** Module renamed to `Feynfeld`;
  package entry moved from `src/v2/FeynfeldX.jl` to `src/Feynfeld.jl` with `v2/`-prefixed
  includes; 71 source/test/script files updated; qgraf submodule's `Main.FeynfeldX.X`
  references collapsed into a proper `import ..Feynfeld: ...` import list (15 new names
  added: `DimD`, `LorentzIndex`, `Spinor`, `feynman_rules`, `propagator_num`,
  `vertex_factor`, `u`/`v`/`ubar`/`vbar`, `spin_sum_amplitude_squared`,
  `spin_sum_interference`, `_degree_partitions`, `_expand_external_fields`,
  `_expand_model_for_diagen`). Tests use `using Feynfeld` (no manual `include`).
- **Move 1.3 (`feynfeld-8cgv`, test orchestration).** `test/v2/runtests.jl` now wires:
  20 original + 5 missing core (`test_diagram_gen`, `test_vertex_arity`,
  `test_qcd_4gluon`, `test_qcd_ghost`, `test_ee_ww_grozin`) + loop-include all 30
  `test/v2/qgraf/test_*.jl` + loop-include `test/v2/munit/test_*.jl` (was run separately
  before). Added thin `test/runtests.jl` top-level forwarder so `Pkg.test()` works.
- **Test count:** 605 → **1327 (1323 pass + 4 `@test_broken` pre-existing dedup bugs)**,
  single-process run ~6 min. `Pkg.test()` verified.
- **Planning info preserved.** Created bd memory `session-28-master-plan` (Move 1/2/3
  structure recovered from the crashed JSONL via Sonnet subagent), epic `feynfeld-8jmm`
  (Move 1 umbrella), epic `feynfeld-0e1t` (Move 3 deferred — "Tobias: I don't care about
  public API yet"). Bead `feynfeld-4206` notes updated with the Move 1.4a-e ambiguity
  (Session 8 closures stand unless user says otherwise).

### Session 29 remaining / out of scope

- **Move 1.4a-e** (5 factory refixes) — NOT executed. Original beads (60n/6mf/blz/3b3/023)
  closed Session 8 with concrete resolution notes; current code's factories still return
  Unions but that may be by-design trade-off. See bead `4206` notes before touching.
- **Move 1.4f** (@inferred regression test file) — still open.
- **Move 1.6** tutorial (`lj1`) — still open.
- **Move 2** (Phase 18b completion) — unchanged. Start at `feynfeld-vjw9`.
- **Move 3** (Process abstraction) — deferred (epic `0e1t`, P4).

### Session 29 addendum: vjw9 attempt + findings (bead still open)

Spent the latter part of Session 29 on `feynfeld-vjw9` (Phase 18b-1a Bhabha
orbit-rep dedup). Red-green TDD with physics ground truth from Session 22
Phase 17a VERDICT and `test/v2/test_bhabha.jl` handbuilt reference.

- **RED test landed**: `test/v2/qgraf/test_phase18b1_multi_orbit.jl` — 2 pass
  (qgraf reference count_diagrams_qg21=2, ee→μμ regression single-orbit) +
  2 `@test_broken` (Bhabha `n_emissions==2`; |M|² == T_tt + T_ss − 2·T_int).
  Handbuilt reference cites `FeynCalc Examples/QED/Tree/Mathematica/ElAel-ElAel.m:98-99`.

- **Three fix attempts, all blocked before green**:
  1. **Signature-based dedup** (`emission_orbit_signature`) — expected orbit-
     invariance via lex-min over autos of `(auto·ps1, _diagram_sig(auto, pmap))`.
     Failed empirically: ee→μμ (1 orbit Burnside) split into 2; Bhabha (2)
     split into 5. The theoretical cocycle identity between `_diagram_sig` and
     a forward pmap action doesn't hold through sorting.
  2. **Direct orbit membership** (`same_emission_orbit`, kept in `audition.jl`,
     exported from `QgrafPort`) — correctly identifies ee→μμ's 8 emissions as
     1 orbit (matches Burnside). Bhabha's 16 emissions split into 4 instead
     of 2; some pmap-action identifications still missing, likely around
     slot-permutation-within-vertex canonicalisation.
  3. **Burnside weighting all emissions × |Stab|/|G|** — paper arithmetic
     verifies |M_t|² + |M_s|² − 2·T_int for Bhabha. **Actual run errors** at
     `interference.jl:46 _find_line_by_bar_mom`: "No line with bar momentum p2"
     on cross-bundle terms. Exactly the secondary bug Session 27 cited when
     picking Strategy B over A ("inter-bundle momentum-label matching inside
     spin_sum_interference ... breaks for automorphic relabelings").

- **Pipeline reverted** to Session 27 `is_emission_canonical` filter so
  Phase 18a regression stays green (ee→μμ pipeline ≡ handbuilt test passes).
  `same_emission_orbit` utility kept in audition.jl + QgrafPort exports —
  it's the useful half of this session's work.

- **Two clear next-attempt paths** (both unlock the `@test_broken` tests):
  - **Option A/C** (bead's naming): finish `same_emission_orbit` for Bhabha.
    Inspect why 4 instead of 2 orbits — likely needs canonical slot-permutation
    within each vertex when comparing the pmap multisets. ~10-30 LOC more.
  - **Option B**: make `spin_sum_interference._find_line_by_bar_mom` robust
    to auto-relabelings — canonicalise bar_mom per bundle before the match.
    Unlocks the already-working Burnside-all arithmetic. ~80-120 LOC
    (Session 27 estimate; tracked in bead `feynfeld-rj1l`).

  Current test suite: 1327 pass (was 1327 pre-attempt) + 2 new `@test_broken`.
  No regressions.

## OPEN FOLLOW-UPS SURFACED SESSION 29

- **`grind/run_v2_tests.sh` is now memory-risky** under concurrent WSL load.
  It pre-dates `test/v2/runtests.jl`. Options: (a) deprecate and document
  `runtests.jl` as canonical, (b) rewrite grind script to reuse one Julia
  process (`julia -e 'include("test/v2/runtests.jl")'`), (c) leave as-is
  and add a banner. No bead filed yet — Tobias to decide.
- **`Pkg.test()` still points at the deleted `test/runtests.jl`.** The
  package's `test/runtests.jl` entry is gone; `Pkg.test()` will now fail.
  Needs either a new top-level `test/runtests.jl` that forwards to
  `test/v2/runtests.jl`, or re-point via `Project.toml`. Part of bead
  `feynfeld-qyu` (FeynfeldX → Feynfeld rename) — the two should land together.
- **qgraf tests (30 files in `test/v2/qgraf/`) still not wired into
  `test/v2/runtests.jl`.** Pre-existing gap, noted in Session 28 stocktake.
- **Previously planned deletions NOT executed this session** (were in the
  original tzgc scope but low-priority):
  - `JULIA_PATTERNS.md` (verbatim duplicate of CLAUDE.md §6) — kept for now.

## SESSION 29 HANDOFF — what the next agent should do

1. **Default next work**: `feynfeld-vjw9` (Phase 18b-1a orbit-rep dedup,
   blocker for Bhabha acceptance). See Session 27 block below for full context.
   File: `src/v2/qgraf/audition.jl:69-95`.
2. **Or tidy the deletion tail**: close `feynfeld-qyu` (FeynfeldX → Feynfeld
   rename) and restore a working `Pkg.test()` entry point. Single-session job.
3. **Or beads hygiene**: walk the 67 open + 6 in-progress beads and decide
   keep/defer/close using the stocktake.
4. **Aggregate test command** (use this, not the grind script):
   `julia --project=. test/v2/runtests.jl`
5. **WSL OOM safety**: if the system feels heavy, check for other WSL Julia
   processes before running large test loops. Single-process runtests.jl is
   the safe default.
6. **Session close protocol**: `bd dolt push` + `git push`.

---

# HANDOFF — 2026-04-17 (Session 28: Full-repo stocktake)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE (Session 28 updates)

1. Read `CLAUDE.md` first. Then this Session 28 block. Then if working code, read
   `reviews/stocktake_2026-04-17/` — six summaries covering every file in the repo:
   - `01_algebra.md` — Layer 4 (20 files, ~2,040 LOC)
   - `02_model_rules_diagrams.md` — Layers 1-3 (15 files, ~1,307 LOC)
   - `03_integrals_evaluate.md` — Layers 5-6 (18 files, ~1,403 LOC)
   - `04_qgraf_port.md` — qgraf port (15 files, 3,583 LOC)
   - `05_tests.md` — 62 test files, 301 @test, 5 @test_broken
   - `06_periphery.md` — v1 frozen, scripts, grind, reviews, refs
2. Current blocker unchanged: `feynfeld-vjw9` (Phase 18b-1a orbit-rep dedup,
   `audition.jl:69` rejects one Bhabha orbit). Next agent starts there unless
   Tobias redirects. See Session 27 block below for full debugging context.
3. No code shipped Session 28. Scope was strictly read-only survey +
   documentation. `bd ready` / `bd list --status=in_progress` unchanged.

## SESSION 28 TIMELINE — full-repo stocktake, no code

1. Tobias: "time to do a stocktake … read the *entire* codebase, generate
   complete documentation of current state … then select files to look at more
   closely … once you truly understand the project report back". Goal stated:
   afterwards we reevaluate which beads to retain.
2. Orchestrated 6 parallel read-only Explore agents (per memory
   `feedback_parallel_agents`: Rule 9 is Julia-only; research agents parallel-OK).
   Each wrote one summary `.md` to `reviews/stocktake_2026-04-17/`. Central-
   summaries pattern to avoid flooding the main context.
3. One agent (periphery) couldn't write because Explore is read-only;
   re-dispatched as general-purpose agent with the same scope. All six landed.
4. Agent 3 (integrals/evaluate) wrote to a top-level mis-path;
   moved into the stocktake dir.
5. Read all six summaries myself to build a complete picture.
6. Drilled directly into 9 high-signal files to verify the agent summaries:
   `src/v2/FeynfeldX.jl` (module entry, 158 LOC, 53 includes, 120+ exports),
   `src/v2/cross_section.jl` (solve_tree + solve_tree_pipeline, 196 LOC),
   `src/v2/qgraf/audition.jl` (the vjw9 blocker file — is_emission_canonical at
   line 69-86), `src/v2/qgraf/burnside_combine.jl` (Session 27's new file, 82 LOC),
   `Project.toml` (5 deps: Combinatorics / LinearAlgebra / PolyLog / QuadGK /
   TensorGR), `src/v2/DESIGN.md`, `src/v2/VERTICAL_PLAN.md`, `SPIRAL_9_PLAN.md`,
   `test/v2/runtests.jl` (20/25 core tests orchestrated).
7. Queried beads state directly (not via summaries, for current truth):
   264 total, 67 open, 7 in-progress, 12 blocked, 190 closed.
8. Reported back; Tobias approved and asked for handoff + commit + push.
9. Wrote this Session 28 block. Added allow-list entry for
   `reviews/stocktake_*/` in `.gitignore` so stocktake snapshots survive
   future fresh clones (pattern changed from `reviews/` to `reviews/*` +
   `!reviews/stocktake_*/` — the bare `reviews/` form blocked re-inclusion
   per git's gitignore rules).

## SESSION 28 ACCOMPLISHMENTS — repo-wide understanding, zero code

- **Every file in the repo accounted for** across the six stocktake summaries.
  v2 source: 53 `.jl` in `src/v2/` (20 algebra + 15 model/rules/diagrams +
  18 integrals/evaluate) + 15 `.jl` in `src/v2/qgraf/` + `FeynfeldX.jl` +
  two in-tree design docs (DESIGN.md, VERTICAL_PLAN.md). v2 tests: 27 main
  + 5 munit + 30 qgraf + runtests.jl. v1 frozen: 28 src + 21 test. Periphery:
  4 scripts, 15+ grind files, 6 architecture reviews + 14 research/ids notes,
  7 refs/ subdirs, 21 papers, 7 top-level docs.
- **Stocktake directory committed** (~70 KB of summaries) so future agents can
  pick up these learnings after a fresh clone instead of re-running six agents.
- **Beads decision layer ready**: with this overview in hand we can walk the
  67 open + 7 in-progress beads category by category and decide keep/defer/close.

## STOCKTAKE FINDINGS — the shape of the repo

### Size and language

**Active (v2):** ~10,300 LOC source across 69 `.jl` files in one module
(`FeynfeldX`); ~7,600 LOC tests across 62 files with 301 @test assertions and
5 @test_broken. 5 Julia deps (Combinatorics, LinearAlgebra, PolyLog, QuadGK,
TensorGR — the last is declared but unused). Package name is still `Feynfeld`
but the module loaded is `FeynfeldX` (bead feynfeld-qyu tracks the rename).

**Frozen (v1):** 28 source files ~3,005 LOC + 21 test files ~2,344 LOC in
`src/algebra/`, `src/integrals/`, `test/algebra/`, `test/integrals/`. Every v1
file has a v2 counterpart except 6 deliberate non-ports (dirac_equation,
dirac_order, dirac_simplify, dirac_scheme, minkowski TensorGR bridge,
feynamp_denominator). `src/{model,rules,diagrams,evaluate}/` and their test
counterparts are empty — dead module scaffolds. `JULIA_PATTERNS.md` (88 LOC)
is a verbatim duplicate of CLAUDE.md §6 / PRD §6.

### What's strong

- **Layer 4 (Algebra) is excellent.** Parametric `Pair{A,B}`, `DiracGamma{S}`,
  `Spinor{K}`, Dict-based AlgSum, `DimPoly` coefficients, all dispatch-based.
  The core architecture validated by the Session 8 six-agent review and by
  Phase 18a's acceptance test: `solve_tree_pipeline(ee→μμ tree massless)` ≡
  `solve_tree(...)` symbolically. The "coefficient type IS the architecture"
  insight in DESIGN.md stays correct — DimPoly eliminated ~150 LOC of v1 glue.
- **Layer 5 (Integrals) is feature-complete for ee→μμ NLO box.** PaVe types,
  B₀ QuadGK, C₀ COLLIER / C0p0 analytical / QuadGK fallback chain, D₀ COLLIER-
  only (triple-nested closure causes JIT explosion, justified), TID rank 0-2.
  50 PaVe tests + 23 D₀ tests cross-validated against LoopTools.
- **qgraf port is deep and faithful.** Phase 17 dedup bug fixed
  (474 → 465 topologies for φ³ 2L via full-permutation Knuth Alg L per equiv
  class). Golden master 95/104 (9 remaining are filter ports + 2 known FAIL).
  Cleanroom `ALGORITHM.md` + per-function `qgraf-4.0.6.f08:XXXX` citations.
  Grind infrastructure in `grind/` allows direct qgraf ↔ Julia trace diffing.

### What's wobbly

- **Pipeline coverage is thin.** Only ee→μμ tree runs end-to-end through
  Layers 1-6. Compton / Bhabha / qq̄→gg / ee→W+W- are hand-built in test files
  (exactly the PIPELINE PRINCIPLE violation the PRD + SPIRAL_9_PLAN.md flag).
  Phase 18b is the fix — 8 sub-tasks wired under epic `feynfeld-xa7s`.
- **Current blocker unchanged from Session 27:** `feynfeld-vjw9`.
  `is_emission_canonical` at `audition.jl:69` rejects one of Bhabha's two
  orbits — canonical-pmap invariant doesn't hold under qgen's flavor
  assignment. HANDOFF Session 22 Phase 17a VERDICT called this the
  Strategy-C under-count. Blocks Phase 18b-1 Bhabha validation.
- **Type instabilities from Session 8 review all still open.**
  `momentum_sum()`, `gamma_pair()`, `pair()` factories return unions;
  `spin_sum.jl:117` + `expand_sp.jl:41,55,69` have `Tuple{Any,...}`;
  `QEDModel.params::Dict{Symbol,Any}` unused; `_COLOUR_DUMMY_COUNTER` global
  mutable. Listed in `reviews/ids_types.txt` (7 beads).
- **MUnit coverage 9 %** (5/60 FeynCalc functions at ≥5 tests). γ5 algebra,
  Eps contraction, DiracEquation, DiracSimplify missing → blocks chiral/EW.
  Spiral 8 was rescoped Session 24 away from "MUnit mop-up"; the MUnit backlog
  (9 P1/P2 beads: DiracTrick 5 batches, DiracTrace 58 tests, EpsContract
  41 tests, SUNTrace, SUNSimplify, ExpandScalarProduct, PolarizationSum,
  Contract) remains.
- **`test/v2/runtests.jl` is incomplete.** 20/25 core tests orchestrated;
  the 30 qgraf tests aren't wired in. `grind/run_v2_tests.sh` covers the main
  25 via separate julia invocations but still skips qgraf. Missing from
  runtests.jl: test_diagram_gen, test_vertex_arity, test_qcd_4gluon,
  test_qcd_ghost, test_ee_ww_grozin.

### Deferrals explicitly marked in code

Phase 18b — 8 sub-tasks, ~560 LOC estimated (see HANDOFF Session 25 table):
18b-1 Burnside multi-orbit (skeleton landed Session 27, needs vjw9 fix),
18b-2 composite-mom fermion propagators, 18b-3 multi-vertex fermion lines
(Compton tree, fermion loops), 18b-4 4-vertex gggg Lorentz, 18b-5 boson
polarisation, 18b-6 symbolic mass, 18b-7 coupling assignment, 18b-8 validation
(Compton, Bhabha, qq̄→gg, ee→W+W-).  Phase 18c (1-loop bridge) blocked on 18b.
Spiral 10 (`feynfeld-4q5`, ee→μμ NLO box via pipeline) blocked on 18c.

### Beads landscape (264 total)

- **67 open** — 19 P1 (Phase 18b sub-tasks + vjw9 + MUnit DiracTrick batches
  + Eps contraction bug + Spiral 10), 31 P2 (MUnit, architectural cleanups,
  golden-master gaps, performance, 2 epics), 8 P3 (ULDM application epic,
  registry, FF library C₀ port), 3 P4 (world-class diagram gen epic, ghost
  fields, native C₀).
- **7 in-progress** — 5 P1 (Phase 14 filters, Phase 17 pipeline swap,
  Phase 18b-1 Burnside, MUnit DiracTrace 58 tests, EpsContract 41 tests),
  2 P2 (Klein-Nishina, Phase 15 symmetry factor).
- **12 blocked**, **190 closed** — the history of Spirals 0-7, Phases 10-18a.

### Periphery findings (recommendations from 06_periphery.md)

Safe to delete when Tobias wants (~5,400 LOC, zero capability loss):
`src/algebra/` (23 files, 2,930 LOC), `src/integrals/` (4 files, 347 LOC),
`src/Feynfeld.jl` (75 LOC), empty `src/{model,rules,diagrams,evaluate}/`,
`test/algebra/` (18 files) + `test/integrals/` (2 files) +
`test/test_ee_mumu.jl` + `test/runtests.jl`, empty
`test/{model,rules,diagrams,evaluate,references}/`, `JULIA_PATTERNS.md`.
Migration sequence in `reviews/stocktake_2026-04-17/06_periphery.md` §10.
**Open question for next session:** re-point `Pkg.test()` at
`test/v2/runtests.jl` and promote `FeynfeldX` → `Feynfeld` (bead feynfeld-qyu).

## SESSION 28 HANDOFF — what the next agent should do

1. **Read the stocktake** (`reviews/stocktake_2026-04-17/01..06`) + this block
   + the Session 27 block below. You will have a complete repo-wide picture.
2. **Then pick a direction:**
   - **Default (unchanged from Session 27):** start on `feynfeld-vjw9` orbit-
     rep dedup. File is `src/v2/qgraf/audition.jl:69-86`. See Session 27
     "Three resolution options" below (A/B/C). Fastest path to unblocking
     Phase 18b-1.
   - **Alt A (beads hygiene):** Tobias may want to walk the 67 open + 7
     in-progress beads and decide keep / defer / close. The stocktake enables
     this — every bead can be evaluated against current-code reality. Do this
     BEFORE starting new 18b work if Tobias asks.
   - **Alt B (v1 deletion):** execute the migration in §10 of
     `06_periphery.md`. ~5,400 LOC removed, `Pkg.test()` re-pointed at
     `test/v2/runtests.jl`, FeynfeldX → Feynfeld rename (bead feynfeld-qyu).
     Single-session job. Close 6+ stale beads that reference v1.
   - **Alt C (type-instability cleanup):** work through
     `reviews/ids_types.txt` — 7 beads from the 2026-03-29 Session 8 review
     still open. Low-risk, high-leverage; removes the CLAUDE.md "MUST FIX"
     list. Each fix is ~10-20 LOC.
3. **Session close protocol** as always: `bd dolt push` + `git push`.

---

## START HERE (Session 27 updates)

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms.
   **Rule 5 note**: 3-research-agents is for older Claude models; Opus 4.7 does
   research by direct reading. Reviewer agent at the end stays mandatory.
2. Run `bd ready` to see available work. **Top priority this session**: `feynfeld-vjw9`
   (Phase 18b-1a orbit-rep dedup, blocker for Bhabha acceptance).
3. Phase 18a regression still green:
   `julia --project=. test/v2/qgraf/test_phase18a_pipeline.jl` → 1/1 ✓
   `julia --project=. test/v2/qgraf/test_solve_tree_pipeline.jl` → 3/3 ✓
4. Phase 18b-1 **skeleton landed** (Session 27) — `burnside_combine.jl` +
   `solve_tree_pipeline` now uses Burnside combine with canonical filter.
   Works for ee→μμ; fails for Bhabha (see blocker below).
5. **Beads created this session**:
   - Epic `feynfeld-xa7s` = Phase 18b umbrella (8 sub-tasks wired with deps)
   - `feynfeld-ewgw` = 18b-1 (claimed, in_progress, blocked by vjw9)
   - `feynfeld-h3pb`, `feynfeld-a7f2`, `feynfeld-m4o8`, `feynfeld-awtt`,
     `feynfeld-feen`, `feynfeld-5d1k`, `feynfeld-4xrh` = 18b-2..8 (open)
   - `feynfeld-rj1l` = Option B best-in-class `InverseSP` factor (future)
   - `feynfeld-vjw9` = 18b-1a orbit-rep dedup (**NEXT AGENT STARTS HERE**)
6. If continuing Phase 18b: start with `feynfeld-vjw9`, unblock 18b-1, then
   write `test/v2/qgraf/test_phase18b1_multi_orbit.jl` (Bhabha acceptance).

---

## SESSION 27 TIMELINE — Phase 18b kickoff

1. Onboarding: read CLAUDE.md, HANDOFF.md, Feynfeld_PRD.md. Internalised the
   pipeline principle, the 12 rules, Session 26 DECISION POINT (A/B/C/D).
2. Recommended **Option A (Phase 18b — tree deferrals)** per HANDOFF Session 25
   rationale: 18a-9 proved the bridge; 18b lifts the artificial scope
   restrictions to make it useful for the full tree-level Standard Model.
3. Tobias: "proceed as you suggest, phase 18b". Created the beads planning
   layer:
   - Epic `feynfeld-xa7s` (Phase 18b umbrella)
   - 8 sub-tasks `ewgw` (18b-1 Burnside), `h3pb` (18b-2 composite-mom fermion
     prop), `a7f2` (18b-3 multi-vertex fermion line, deps on h3pb),
     `m4o8` (18b-4 boson polarisation), `awtt` (18b-5 4-vertex gggg),
     `feen` (18b-6 symbolic mass), `5d1k` (18b-7 coupling assignment),
     `4xrh` (18b-8 validation, deps on all).
4. Rule-5 check: Tobias clarified **"3 agent rule is for older models of
   Claude"**. Saved to memory (feedback_core_rules_discipline.md updated).
   Opus 4.7: read source directly + reviewer agent at end.
5. Read handbuilt path to scope 18b-1: `interference.jl` (spin_sum_interference,
   _cross_line_trace), `spin_sum.jl` (spin_sum_amplitude_squared, _single_line_trace,
   _conjugate_gammas with :mu → :mu_ relabel), `amplitude.jl`
   (:mu_<channel> naming, _fermion_line_chain), `test/v2/test_bhabha.jl`
   (|M|² = (1/4)(T_tt/t² + T_ss/s² − 2·T_int/(s·t)) handbuilt).
6. **Open question surfaced**: AlgSum has no inverse-denom factor — cannot
   symbolically represent 1/pair(q,q). Proposed **Option A** (return trace-only
   AlgSum; caller applies 1/denom) vs **Option B** (introduce `InverseSP` factor).
   Tobias: "sounds good, but I do want best in class solutions eventually, so
   make sure option B is recorded as a bead". Created `feynfeld-rj1l`.
7. **Implemented Phase 18b-1 skeleton** (Option A scope):
   - New file `src/v2/qgraf/burnside_combine.jl` (~80 LOC) —
     `combine_m_squared_burnside(bundles, weights)` + `_pair_trace` helper.
     Uses `spin_sum_amplitude_squared` for diagonal, `spin_sum_interference`
     for off-diagonal, Burnside weights w_i·w_j, and fermion signs
     bi.fermion_sign · bj.fermion_sign.
   - `src/v2/qgraf/QgrafPort.jl` — include + export `combine_m_squared_burnside`.
   - `src/v2/cross_section.jl::solve_tree_pipeline` — replaced the Phase 18a
     `bundles[1]` single-orbit shortcut with the Burnside combine. Added
     canonical-rep filter (`is_emission_canonical`) to the emission loop so
     weights collapse to 1 per orbit (see Session 27 blocker below).
     Return tuple extended: `(amplitude_squared, n_emissions, orbit_denoms)`.
     196 LOC total — under Rule 11 ~200 ceiling.
8. **Phase 18a regression verified green** after wiring:
   `test_phase18a_pipeline.jl` 1/1 ✓; `test_solve_tree_pipeline.jl` 3/3 ✓.
9. **Bhabha blocker surfaced** — Session 27 key finding:
   `solve_tree_pipeline(Bhabha ee→ee tree)` reports `n_emissions=1`, but
   `count_diagrams_qg21(qed_model, [:e,:e_bar], [:e,:e_bar])` correctly
   returns 2 (s + t orbits). The canonical filter rejects one of the two
   orbits — exactly the Strategy-C under-count bug documented in HANDOFF
   Session 22 Phase 17a VERDICT: "the canonical orbit-rep may be INVALID
   for qgen, so the orbit yields 0 emissions instead of 1 (under-count)."
   Filed as blocker `feynfeld-vjw9` (Phase 18b-1a).
10. Updated `feynfeld-ewgw` (18b-1) notes with current status and dep on vjw9.
11. Stopped here at Tobias's request ("stop at the next most convenient place").
    Machinery in place; Phase 18a still green; Bhabha validation blocked on
    orbit-rep dedup.

## SESSION 27 ACCOMPLISHMENTS

- Phase 18b planning: 9 beads (epic + 8 sub-tasks + Option B retrofit), full
  dependency graph wired (`bd dep add` chain). +1 blocker bead post-debug.
- Phase 18b-1 implementation skeleton: ~80 LOC new, ~10 LOC delta in
  cross_section.jl, QgrafPort exports extended.
- Phase 18a regression stays green (ee→μμ canonical-filter-compatible).
- Memory update: `feedback_core_rules_discipline.md` now reflects Opus 4.7
  relaxation of the 3-agent rule.

## NEXT AGENT: START WITH `feynfeld-vjw9` (Phase 18b-1a)

### What's broken and why

`solve_tree_pipeline(Bhabha)` returns `n_emissions=1`. Bhabha has 2 orbits
(s-channel photon annihilation + t-channel photon exchange). My canonical
filter via `is_emission_canonical` (audition.jl:69-95) discards one of them.
This is the **Strategy-C under-count bug** documented in
HANDOFF Session 22 Phase 17a VERDICT: Strategy C "pre-filters ps1 to
orbit-reps assuming the rep is qgen-valid, which it may not be."

Evidence reproducible in ~10s:
```julia
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
using .FeynfeldX.QgrafPort: count_diagrams_qg21
p1=Momentum(:p1); p2=Momentum(:p2); k1=Momentum(:k1); k2=Momentum(:k2)
prob = CrossSectionProblem(
    qed_model(m_e=:zero, m_mu=:zero),
    [ExternalLeg(:e, p1, true,  false), ExternalLeg(:e, p2, true,  true)],
    [ExternalLeg(:e, k1, false, false), ExternalLeg(:e, k2, false, true)],
    10.0,
)
println("count_diagrams_qg21 (Strategy A Burnside) = ",
        count_diagrams_qg21(qed_model(), [:e,:e_bar], [:e,:e_bar]; loops=0))
println("solve_tree_pipeline n_emissions (canonical filter) = ",
        solve_tree_pipeline(prob).n_emissions)
'
```
Expected: 2 vs 1.

### Why I used the canonical filter instead of Burnside-all

Burnside-all (sum every orbit member with weight 1/|Orbit|) is robust for
COUNTING (Strategy A chose this) but breaks multi-orbit amplitude summation:
`spin_sum_interference` (interference.jl:102 `_find_line_by_bar_mom`) keys
matching on bar-momentum names. Within one orbit, different members have
automorphic momentum relabelings (build_externals at vertex_assemble.jl:141
binds field → physical_moms[i] via pmap[i,1], and pmap varies across orbit
members) → loop-close fails. Canonical filter sidestepped this by keeping
one rep per orbit, but hits the Strategy-C validity bug.

### Three resolution options (documented in feynfeld-vjw9)

**(A) Fix is_emission_canonical** to fall back to lex-next-smallest when
lex-smallest is qgen-invalid. Narrow fix to the Strategy-C bug; preserves
my current solve_tree_pipeline wiring. ~20-30 LOC in audition.jl.

**(B) Switch to Burnside-all + canonical relabeling** — sum every orbit
member with weights 1/|Orbit| (removes Strategy-C bug), but canonicalise
bar_mom per bundle before calling `spin_sum_interference`. Requires
writing a per-bundle momentum relabeler. ~80-120 LOC, more general.

**(C) Hybrid** — use Strategy-A counting to enumerate orbits, then use
a signature (hash of canonical pmap) to group and pick one qgen-valid
member per orbit. ~50 LOC, moderate complexity.

**Recommendation**: start with (A) — it's the targeted fix for a known
bug and unblocks 18b-1 fastest. If (A) reveals deeper issues with the
lex ordering, escalate to (C). (B) is worth doing eventually but is
18b-3-adjacent (when internal fermion propagators need relabeling anyway).

### After vjw9 closes

1. Write `test/v2/qgraf/test_phase18b1_multi_orbit.jl` — Bhabha acceptance:
   - `result.n_emissions == 2`
   - `result.amplitude_squared == handbuilt_trace_only` where
     `handbuilt_trace_only = T_tt + T_ss − 2·T_int` (NO denoms, trace only —
     Option A scope per feynfeld-rj1l). Indices: handbuilt uses `:alpha`/`:alpha_`
     for t-channel and `:beta`/`:beta_` for s-channel. Pipeline uses
     `:mu_l_<edge_id>` naming. If symbolic `==` fails due to index-label
     differences, first try contract+expand on both sides (indices are dummy,
     should normalise away). If it STILL fails, unify naming per
     vertex_assemble.jl:62-69 comment.
   - Reference handbuilt derivation in `test/v2/test_bhabha.jl:63-94`.
2. Once green: close feynfeld-ewgw (18b-1), move to next 18b sub-task.
   Recommend `feynfeld-h3pb` (18b-2 composite-mom fermion prop) next since
   it unblocks 18b-3 (multi-vertex fermion line = Compton tree validation).

### Potential traps

- `emission_amplitude.jl:64-68` comment says "negate outgoing" but the loop
  just copies physical → qgraf_ext_moms. Masked by ee→μμ symmetry. Bhabha
  may surface this. Flag if discrepancies appear.
- Boson Lorentz index divergence (`vertex_assemble.jl:62-69`): pipeline's
  `:mu_l_<edge_id>` vs handbuilt's `:mu_<channel>`. Currently masked by
  in-chain contraction. Bhabha cross-terms may or may not tolerate this.
  In-source comment flags the unification task.
- `spin_sum_interference` (interference.jl:44-45) currently errors on
  multi-term DiracExpr (chiral vertices). QED Bhabha has single-term γ^μ
  vertex — fine. EW will need extension later (18b-4 territory).

## FILES TOUCHED THIS SESSION

| Path | Change | LOC |
|------|--------|----:|
| `src/v2/qgraf/burnside_combine.jl` | NEW | 82 |
| `src/v2/qgraf/QgrafPort.jl` | include + export | +2 |
| `src/v2/cross_section.jl` | solve_tree_pipeline Burnside + canonical filter + 3-field return | ~+10/-20 |
| `HANDOFF.md` | this session | (you're reading it) |

No test files added (18b-1 acceptance test deferred until vjw9 unblocks it).
No regressions in existing tests.

## QUICK COMMANDS

```bash
# Phase 18a regression (should stay green):
julia --project=. test/v2/qgraf/test_phase18a_pipeline.jl
julia --project=. test/v2/qgraf/test_solve_tree_pipeline.jl

# Bhabha blocker reproduction:
julia --project=. -e 'include("src/v2/FeynfeldX.jl"); using .FeynfeldX;
using .FeynfeldX.QgrafPort: count_diagrams_qg21;
p1,p2,k1,k2 = Momentum.((:p1,:p2,:k1,:k2));
prob = CrossSectionProblem(qed_model(m_e=:zero,m_mu=:zero),
  [ExternalLeg(:e,p1,true,false), ExternalLeg(:e,p2,true,true)],
  [ExternalLeg(:e,k1,false,false), ExternalLeg(:e,k2,false,true)], 10.0);
println("orbits: ", count_diagrams_qg21(qed_model(),[:e,:e_bar],[:e,:e_bar]));
println("emissions: ", solve_tree_pipeline(prob).n_emissions)'

# Beads status:
bd ready                # pick up feynfeld-vjw9 first
bd show feynfeld-vjw9   # read the blocker
bd show feynfeld-ewgw   # read 18b-1 current state

# Pre-existing green suites (before touching anything):
julia --project=. test/v2/test_diagram_gen.jl      # 32/32 ✓
./grind/run_v2_tests.sh                            # 26/26 files
```

---

## SESSION 26 TIMELINE (exploration only — no code shipped)

1. Read `HANDOFF.md`, `Feynfeld_PRD.md`, `CLAUDE.md`. Internalised THE
   PIPELINE PRINCIPLE, the 12 rules, the spiral methodology, and the
   Session 25 NEXT SESSION DECISION POINT (Phase 18b vs 18c vs filter ports).
2. Tobias: "I just installed Mathematica, and hence wolframscript. Can you
   gently explore the capabilities and whether you can use wolframscript
   and FeynCalc etc. to generate golden masters for other parts of the
   pipeline." Shifted from the pending Phase 18b to an exploratory scoping
   session (no Julia code touched — safe read-only exercise).
3. Environment verification (all green):
   - `/usr/bin/wolframscript`, WolframScript 1.13.0 for Linux x86 (64-bit)
   - Mathematica 14.3.0 for Linux x86 (July 2025 build)
   - Install: `/usr/local/Wolfram/Wolfram/14.3`; user base: `~/.Wolfram`
   - FeynCalc source tree: `refs/FeynCalc/FeynCalc/` (entry `FeynCalc.m`)
4. Load-pattern established (**no paclet install required**):
   ```bash
   wolframscript -code 'PrependTo[$Path, "/home/tobiasosborne/Projects/Feynfeld.jl/refs/FeynCalc"]; Needs["FeynCalc`"]; ...'
   ```
   Confirmed: FeynCalc 10.2.0 (dev version, 22cc5e08, 2026-03-25) loads in
   ~5 s. The `$FCTraceNames` multi-context shadow warning is cosmetic.
5. Capability probe #1 — Dirac traces:
   - `DiracTrace[GA[μ].GA[ν]]` → `4 Pair[LorentzIndex[μ], LorentzIndex[ν]]`
   - `DiracTrace[GA[μ].GA[ν].GA[ρ].GA[σ]]` → canonical 3-term g·g·g form
   - `DiracTrace[GA[μ].GA[ν].GA[ρ].GA[σ].GA[5]]` → `-4I Eps[LorentzIndex[...]]`
6. Capability probe #2 — PaVe reduction (Layer 5 hole):
   - `PaVeReduce[B1[p², m0², m1²]]` → full Denner Fortschr. Phys. 41 (1993)
     Eq. (4.18) analytically, in `PaVe[0, ..., {...}]` form.
   - `Tdec[{{q, μ}}, {p}]` → tensor decomposition basis (symbolic
     Passarino-Veltman invariants).
7. Scoped 5 golden-master surfaces ranked by leverage:
   1. **PaVe reduction library** (Layer 5) — B/C/D ij, D ijk — direct oracle
      for the PV reduction layer. Highest physics ROI; narrow symbolic
      surface. **Recommended first surface.**
   2. **γ5 / Eps golden suite** — unblocks Spiral 8 (the pending spiral).
      HVBM vs NDR selectable via `$BreitMaison`.
   3. **Dirac trace battery** — up to 10-gamma traces, direct oracle for v2
      `dirac_trace`; comparable to FeynCalc `Tests/Dirac/DiracTrace.test`.
   4. **End-to-end |M|²** — FeynArts + FCFAConvert + FeynCalc for
      ee→μμ, Compton, Bhabha, qq̄→gg (cross-validates Phase 18a bridge).
   5. **SUN color golden** — `SUNSimplify` traces for QCD (Spiral 3+).
8. Proposed workflow (documented, not implemented):
   ```
   scripts/golden_master_<surface>.wls   (wolframscript generator)
       ↓ emits
   test/v2/golden/test_<surface>.jl      (Julia @testset, citations)
   ```
   Each `.wls` emits FeynCalc `InputForm` + hand-mapped Feynfeld expression,
   tagged with FeynCalc MUnit test IDs where applicable. Both the generator
   and its output are checked in so results are reproducible.
   Estimated per-surface cost: ~150 LOC Wolfram + ~100 LOC Julia translator
   + ~40 golden tests. ~50 LOC of symbolic-normalisation per surface for
   the `Pair[LorentzIndex[μ],...]` ↔ `pair(μ,ν)` mapping.
9. Tobias: "yea give it a go". Kicked off parallel research per Rule 5:
   - Read `src/v2/pave.jl` (v2 PaVe type: `PaVe{N}` parametric, named
     constructors A0/B0/B1/B00/B11/C0/C1/C2/D0/D1/D2/D3, canonical field
     ordering `sort(indices)`).
   - Attempted `ls refs/FeynCalc/FeynCalc/Tests/LoopIntegrals/` → **path
     wrong**, no `Tests/` directory exists under the FeynCalc paclet.
     The real FeynCalc MUnit tests live at top level `refs/FeynCalc/Tests/`
     (parallel to `FeynCalc/`); NOT yet confirmed this session.
   - Listed `test/v2/` — 27 test files + `munit/` + `qgraf/` + `runtests.jl`.
     No `golden/` directory yet; existing `test_pave.jl` is the neighbour.
10. **Interrupted before** creating the beads epic or drafting the
    wolframscript generator or the translator. Session ended here at
    Tobias's request. No code changes; only `.beads/backup/backup_state.json`
    was touched by `bd remember` at session end.

## SESSION 26 ACCOMPLISHMENTS — scoping complete, implementation not started

- **Environment verified**: wolframscript 1.13 + Mathematica 14.3 + FeynCalc
  10.2.0 load via `PrependTo[$Path, ...] + Needs["FeynCalc\`"]` with no
  install step. 5 capability probes all green.
- **Strategic landscape mapped**: 5 golden-master surfaces ranked by ROI;
  PaVe reduction (Layer 5) chosen as recommended first surface on leverage
  and narrow symbolic surface grounds.
- **Workflow template drafted** (documented in this HANDOFF, not yet
  committed as code): `scripts/golden_master_<surface>.wls` → `test/v2/golden/`.
- **Persistent memory saved** via `bd remember --key feyncalc-wolframscript-setup`
  so the next session loads the verified probe outputs and the load pattern.

## WHAT SESSION 26 ENABLES (when the generator is built)

- **Layer 5 closes faster**: v2 currently has only A0, B0, B1, C0, C1, C2
  hand-coded (plus B00/B11/D tensor stubs per pave.jl:34-66). FeynCalc can
  emit the full Passarino-Veltman set (B_ij, C_ij, C_ijk, D_ij, D_ijk, D_ijkl)
  with textbook-matching analytical reduction. Cuts Phase 18c (1-loop bridge)
  prep work significantly.
- **γ5 convention locked**: Spiral 8 needs HVBM vs NDR choice made and
  frozen. FeynCalc exposes `$BreitMaison` as a single-flag switch; a golden
  suite in both schemes documents the choice unambiguously.
- **Cross-validation target**: Phase 18a ≡ handbuilt is one test. Phase
  18a ≡ FeynCalc FCFAConvert output is a second, independent witness —
  much stronger evidence the bridge is right.
- **Anti-hallucination hardening**: MUnit tests in `refs/FeynCalc/Tests/`
  become mechanically reproducible via the generator, removing the manual
  translation step for the routine-permutation tier (CLAUDE.md §MUnit
  translation protocol 4b).

## NEXT SESSION CONCRETE STEPS (if continuing Option D)

1. **Scout** the actual FeynCalc `Tests/` layout (this session's path guess
   was wrong; the tests almost certainly live at `refs/FeynCalc/Tests/` with
   subdirs `Dirac/`, `LoopIntegrals/`, `Lorentz/`, `SUN/`, …).
2. **Create beads epic** `feynfeld-PHASE19-feyncalc-goldens` with sub-tasks:
   - 19-1: Scout Tests/ layout + pick a first handful of PaVe tests to port.
   - 19-2: Draft `scripts/golden_master_pave.wls` — emits the B1, B00, B11,
     C_ij, D_ij reductions in `InputForm`.
   - 19-3: Draft Julia translator (`PaVe[0, {p²}, {m0², m1²}]` →
     `B0(p2, m02, m12)`; `PaVe[1, ...]` → `B1(...)`; etc.). Round-trip
     test: translator(generator output) should equal Feynfeld's pave.jl
     constructors.
   - 19-4: Emit `test/v2/golden/test_pave_reduction.jl` with first ~20 goldens.
   - 19-5: Wire into `./grind/run_v2_tests.sh`.
   - 19-6: Per Rule 6, rigorous reviewer agent on the generator + translator.
3. **Do not start 19-2 without Rule 5 tiered research** — this is core
   infrastructure (crosses two languages, affects how every future Layer 4/5
   test is validated). Tiered workflow: 3 research + 1 review.
4. Estimated session size: 1 full session for PaVe surface (~300 LOC
   Wolfram + Julia + ~40 goldens), not counting review iteration.

## KNOWN GAPS / RISKS for Option D

- **Symbolic normalisation is per-surface work**. FeynCalc emits
  `Pair[LorentzIndex[μ], LorentzIndex[ν]]`, `SPD[p,q]`, `FVD[p,μ]`,
  `DiracGamma[LorentzIndex[μ]]`, `PaVe[i, {invs}, {masses}]`. Feynfeld
  uses `pair(μ,ν)`, `SP(p,q)`, own `DiracChain`, `PaVe{N}` parametric.
  The translator is not a one-liner; expect ~50-100 LOC per surface.
- **γ5 convention**: FeynCalc defaults to `$BreitMaison=False` (NDR).
  Feynfeld has not yet committed to a scheme. **Decision required**
  before γ5 golden suite lands. Not blocking for PaVe (first surface).
- **Rule 1 still rules**: even with a FeynCalc golden, the test must cite
  the textbook equation that validates it. FeynCalc is an oracle, not a
  primary source. Routine-permutation tier (§MUnit protocol 4b) is ok
  with FeynCalc-only citation.
- **No `JuliaForm`**. Wolfram has `CForm`, `FortranForm`, `TeXForm`, but
  no Julia emitter. Output is via `InputForm` + Julia-side parser, or
  structural walk of the Mathematica expression. The latter is cleaner
  for non-trivial trees.

## SESSION 25 TIMELINE

1. Read `HANDOFF.md`, `Feynfeld_PRD.md`, `CLAUDE.md`. Internalised the
   THE PIPELINE PRINCIPLE, the 12 rules, the spiral methodology, and
   the Session 24 NEXT SESSION DECISION POINT (Options A/B/C).
2. Picked **Option A — Phase 18a tree-level Diagram → AlgSum bridge**
   per the recommendation. The qg21 port was a diagram counter; this
   phase makes it produce evaluable amplitudes.
3. Drafted the granular plan: 10 sub-tasks (18a-1 through 18a-10),
   created beads epic `feynfeld-otgb` + 10 task issues with explicit
   dependency graph (`bd dep add` chain).
4. Per Tobias Rule 5 (core code = 3 research + 1 review), spawned 3
   parallel read-only research agents on qgraf f08:13400-13580 leaf-peel
   (algorithm details), Nogueira 1993 paper (cross-check), and existing
   Julia momentum API surface. Verified all three reports against my
   own direct read of f08:13313-13491.
5. **Phase 18a-1** (`f4563bb`): RED-GREEN TDD test by test. Built
   `route_momenta(state, labels, ext_moms; loop_moms)` returning
   `EdgeMomenta` (per-edge `MomentumSum` + edge-type tag). Tests
   walked: φ→φφ tree, ee→μμ s-channel (p1+p2), ee→μμ t-channel
   (p1+p3), φ³ 1L bubble (chord head-match flip), φ³ tadpole
   (snb edge type), qg21_enumerate integration. Side-fix:
   `MomentumSum == /hash` (default Julia struct == falls back to
   === for Vector-bearing types). Commit ~220 LOC source + 25 tests.
6. **Phase 18a-2** (`85a325d`): `compute_amap(state, labels)` —
   half-edge labelling matrix. External back-write + internal triple-
   case (single edge / self-loop integer division / parallel edge
   backward scan), per qgraf f08:12133-12158 + 12342-12344. RED-GREEN
   tests: φ→φφ, ee→μμ, φ³ bubble parallel edges, φ³ tadpole
   self-loop, comprehensive pairing-invariant battery (4 topologies).
   Side-fix: qgen.jl `vmap`/`lmap` allocation `n×n → n×MAX_V` (latent
   self-loop bug — vdeg can exceed n, surfaced by tadpole test).
   ~80 LOC source + 57 tests.
7. **Phase 18a-3** (`60acb42`): `build_propagators` — per-edge
   propagator factors (Boson `alg(1)` / Fermion `DiracExpr(p̸+m)` /
   Scalar `alg(1)`). Denominator = `pair(mom, mom) − m²`. Tests:
   ee→μμ photon, φ³ scalar, φ³ bubble (2 parallel propagators),
   tadpole self-loop. ~120 LOC + 24 tests.
8. **Phase 18a-4** (`149ba8a`): `build_vertices` — per-vertex
   Lorentz factors. Boson edge index naming `:mu_l_<edge_id>` shared
   between endpoints (Einstein summation auto-contracts at chain
   product). Field canonicalisation strips `_bar` for model dict
   lookup. Tests: ee→μμ γ^μ at both vertices with shared index, φ³
   scalar (no Lorentz). Side-fix: `DiracChain` and `DiracExpr` ==
   /hash methods (same Vector-equality root cause as 18a-1).
   ~120 LOC + 7 tests.
9. **Phase 18a-5** (`1c64820`): `build_externals` — per-external
   spinor/polarisation. Mirrors amplitude.jl `_spinor_and_position`
   dispatch (u/v/ubar/vbar from in/out + antiparticle flags). Boson
   externals deferred (returns nothing). ~80 LOC + 19 tests.
10. **Phase 18a-6** (`70739a2`): `walk_fermion_lines` — pairs each
    internal vertex's 2 fermion half-edges by bar/plain end via
    `build_externals`' position metadata. Tree-only: errors with a
    Phase-18b deferral message if a fermion slot connects to another
    internal vertex (Compton-style internal fermion propagator).
    ~85 LOC + 10 tests.
11. **Phase 18a-7** (`eebf79f`): `emission_to_amplitude` — master
    assembler. Composes 18a-1..6 into `AmplitudeBundle(line_chains,
    amplitude, denoms, fermion_sign, sym_factor, coupling)`. Per-line
    chain construction mirrors amplitude.jl `_fermion_line_chain`.
    Tests: ee→μμ s-channel bundle structure, φ³ scalar bundle (no
    fermion lines → `DiracExpr(alg(1))`). ~135 LOC + 9 tests.
12. **Phase 18a-8** (`40dc142`): `solve_tree_pipeline` — drives the
    qg21 `_foreach_emission` stream into `emission_to_amplitude`,
    picks the first emission's bundle (single-orbit Phase-18a
    shortcut), runs `spin_sum_amplitude_squared` → `contract` →
    `expand_scalar_product`. ~55 LOC + 3 smoke tests.
13. **Phase 18a-9** (`cdb0262`): THE acceptance test —
    `solve_tree_pipeline(qed_model, ee→μμ massless).amplitude_squared
    == solve_tree(...).amplitude_squared` symbolically. PASS first
    try. The pipeline produces the same |M|² as the hand-built path
    after spin-sum/contract/expand_sp. **Phase 18a milestone
    achieved.**
14. **Phase 18a-10** (`cd39db4`): Spawned read-only reviewer agent
    on commits `f4563bb..cdb0262`. Verdict: SHIP-READY with 3
    caveats: (a) momentum.jl LOC limit (297 → split spanning_tree.jl
    off, now 230); (b) boson Lorentz index naming divergence
    (`:mu_l_<edge_id>` vs `:mu_<channel>` — masked by DiracChain
    contraction but flagged for Phase 18b unification with in-source
    comment at vertex_assemble.jl:65); (c) side-fix commits should
    have been standalone (decided to leave bundled — explicit in
    commit messages). HANDOFF updated, beads closed (epic feynfeld-otgb
    + 10 tasks), `git push` + `bd dolt push` clean.

## SESSION 25 ACCOMPLISHMENTS — Phase 18a CLOSED

The qg21 port is no longer just a diagram counter — it now produces
evaluable AlgSum amplitudes. End milestone proven:
**`solve_tree_pipeline(qed_model, ee→μμ massless) ==
solve_tree(qed_model, ee→μμ massless)`** symbolically, after spin-sum
/ contract / expand_sp.

### Phase-by-phase (10 commits f4563bb..cdb0262)

| Phase | What | LOC | Tests | Commit |
|-------|------|----:|------:|--------|
| 18a-1 | leaf-peel `route_momenta` (qgraf f08:13400-13559) | ~220 | 25 | `f4563bb` |
| 18a-2 | half-edge `compute_amap` (f08:12133-12158) | ~80 | 57 | `85a325d` |
| 18a-3 | per-edge `build_propagators` (Boson/Fermion/Scalar) | ~120 | 24 | `60acb42` |
| 18a-4 | per-vertex `build_vertices` (γ^μ + index sharing) | ~120 | 7 | `149ba8a` |
| 18a-5 | per-external `build_externals` (u/v/ubar/vbar) | ~80 | 19 | `1c64820` |
| 18a-6 | fermion-line `walk_fermion_lines` (tree-only) | ~85 | 10 | `70739a2` |
| 18a-7 | master `emission_to_amplitude` → AmplitudeBundle | ~135 | 9 | `eebf79f` |
| 18a-8 | `solve_tree_pipeline` (cross_section.jl wiring) | ~55 | 3 | `40dc142` |
| 18a-9 | symbolic equality test ee→μμ pipeline ≡ handbuilt | ~30 | 1 | `cdb0262` |

Net: ~925 LOC source, 9 new test files, 155 new tests, 1 critical
acceptance test passing.

### Side fixes (bundled in phase commits, all triggered by Phase 18a)

- `MomentumSum == / hash` (types.jl) — default Julia `==` was `===`
  for Vector-bearing struct; added in 18a-1.
- `vmap/lmap` allocation `n×n → n×MAX_V` (qgen.jl) — bug surfaced by
  the φ³ tadpole test (vdeg can exceed n with self-loops); 18a-2.
- `DiracChain == / hash` (dirac.jl) — same Vector-equality issue; 18a-4.
- `DiracExpr == / hash` (dirac_expr.jl) — same; 18a-4.

### New module structure (src/v2/qgraf/)

| File | Purpose |
|------|---------|
| `momentum.jl` (extended) | Phase 16 spanning tree + Phase 18a-1 leaf-peel |
| `halfedge.jl` (new) | compute_amap |
| `propagator_assemble.jl` (new) | build_propagators |
| `vertex_assemble.jl` (new) | build_vertices + build_externals |
| `fermion_line.jl` (new) | walk_fermion_lines (tree-only) |
| `emission_amplitude.jl` (new) | emission_to_amplitude (master assembler) |

### Reviewer findings (post-18a-9)

A read-only review agent flagged 3 caveats on the closed phase. Disposition:

1. **Boson Lorentz index naming divergence** (vertex_assemble.jl:65):
   pipeline uses `:mu_l_<edge_id>`, handbuilt uses `:mu_<channel>`.
   Currently masked by contraction inside DiracChain dot products
   (both produce the same final AlgSum). Will need unification before
   Phase 18b's explicit metric / multi-orbit / polarisation work.
   **Action**: in-source comment added at vertex_assemble.jl:65.
2. **momentum.jl LOC (was 297)**: Split spanning_tree.jl off (69 LOC),
   leaving momentum.jl at 230 LOC — closer to the ~200 rule. Route_momenta
   itself is 150 LOC of dense leaf-peel logic and doesn't split cleanly.
3. **Side-fix commits** (MomentumSum/DiracChain/DiracExpr ==/hash,
   vmap/lmap allocation): each was bundled into the phase commit that
   triggered it, with clear documentation in the commit message. Reviewer
   suggested extracting via rebase. **Disposition**: leave as-is — the
   commit messages are explicit and the changes were all symmetric to
   existing patterns. A future cleanup pass can extract if desired.

### WHAT PHASE 18a ENABLES

#### 1. Architectural — pipeline principle satisfied for tree QED 2→2 boson exchange

Before 18a, every physics process bypassed Layers 1-3 and used
hand-rolled amplitudes (channels.jl, amplitude.jl, build_amplitude
per process). CLAUDE.md's THE PIPELINE PRINCIPLE was aspirational.
After 18a, for the validated subset (ee→μμ tree massless), the full
6-layer pipeline runs end-to-end: Model → Rules → Diagrams (qg21) →
Algebra (AlgSum, DiracExpr) → Integrals (PaVe-ready) → Evaluate
(spin-sum, contract, expand_sp). No bypass.

#### 2. Concretely usable APIs

```julia
# Drop-in replacement for solve_tree (validated symbolic equivalence):
prob = CrossSectionProblem(qed_model(m_e=:zero, m_mu=:zero),
                           [ExternalLeg(:e, p1, true,  false),
                            ExternalLeg(:e, p2, true,  true)],
                           [ExternalLeg(:mu, k1, false, false),
                            ExternalLeg(:mu, k2, false, true)],
                           10.0)
result = solve_tree_pipeline(prob)
result.amplitude_squared isa AlgSum   # spin-summed |M|²
result.n_emissions                    # qg21 emission count
```

```julia
# One-shot "give me the amplitude bundle for this emission":
bundle = emission_to_amplitude(state, labels, ps1, pmap, model;
                                physical_moms, n_inco)
bundle.line_chains   # Vector{DiracExpr}, one per fermion line
bundle.amplitude     # convenience: product of line_chains
bundle.denoms        # Vector{AlgSum} — (p²-m²) per internal propagator
bundle.fermion_sign  # ±1 from qdis_fermion_sign
bundle.sym_factor    # 1/S_local (Rational)
```

The 6 sub-builders (`route_momenta`, `compute_amap`,
`build_propagators`, `build_vertices`, `build_externals`,
`walk_fermion_lines`) are independently usable for diagnostics,
partial assembly, or alternative amplitude conventions.

#### 3. Validation pattern established

`pipeline ≡ handbuilt symbolic AlgSum equality` is now a template
test pattern (`test/v2/qgraf/test_phase18a_pipeline.jl`). Each
hand-built process (Compton, Bhabha, qq̄→gg, vertex_g2, etc.) can
get an analogous test once the relevant 18b deferral is lifted. The
hand-built code becomes a "ground-truth oracle" while the pipeline
catches up — and once oracular tests pass, the hand-built path can
eventually be retired.

#### 4. Downstream pipelines unblocked

- **Layer 5 (PaVe)**: `AmplitudeBundle.denoms` is the list of
  `(p²−m²)` factors a 1-loop variant feeds into Passarino–Veltman
  reduction. Phase 18c (1-loop) becomes a structural extension, not
  a redesign.
- **Layer 6 (cross section, observables)** already consumes AlgSum
  via `evaluate_m_squared`, `dsigma_domega`, `evaluate_numeric` —
  these now work transparently on pipeline output. Tree-level
  cross sections via the pipeline work end-to-end (modulo what's
  deferred to 18b).

#### 5. Agent-facing value (PRD §1.2 endgame)

The PRD vision: "Claude reads Lagrangian → returns σ_NLO". Until
18a, every spoke required Claude to write custom Layer-3 code. Now:
for any process whose deferrals 18b lifts, the agent invokes
`solve_tree_pipeline(CrossSectionProblem(...))` and the rest is
automatic. This is the first session where the bridge to
"agent-driven physics" exists in code, not just documentation.

#### 6. Module surface added (src/v2/qgraf/)

Exported via QgrafPort:
- Types: `EdgeMomenta`, `InternalEdge`, `Propagator`, `ExternalFactor`,
  `FermionLine`, `AmplitudeBundle`
- Functions: `route_momenta`, `compute_amap`, `build_propagators`,
  `build_vertices`, `build_externals`, `walk_fermion_lines`,
  `emission_to_amplitude`, `_foreach_emission` (re-export of audition.jl)

Exported via FeynfeldX (Layer 6):
- `solve_tree_pipeline`

### Phase 18b roadmap (concrete)

To complete tree-level Standard Model coverage:

| Sub-task | What to lift | Where | Est. LOC |
|----------|--------------|-------|---------:|
| 18b-1 | Multi-orbit Burnside summation | cross_section.jl:155 | ~50 |
| 18b-2 | Fermion propagator with composite momentum | propagator_assemble.jl:88 | ~30 |
| 18b-3 | Multi-vertex fermion-line traversal | fermion_line.jl:55 | ~120 |
| 18b-4 | Boson polarisation (external gluons) | vertex_assemble.jl ext branch | ~80 |
| 18b-5 | 4-vertex (gggg) Lorentz factor | vertex_assemble.jl:30 | ~60 |
| 18b-6 | Symbolic mass placeholders | propagator_assemble.jl:75 | ~40 |
| 18b-7 | Coupling assignment (e², g_s², etc.) | emission_amplitude.jl:140 | ~30 |
| 18b-8 | Validation: Compton, Bhabha, qq̄→gg, ee→W+W- | test/v2/qgraf/ | ~150 |

Total estimate: ~560 LOC, 1-2 sessions. Each unlocks one or more
hand-built processes for symbolic-equivalence cross-validation.

### Phase 18c sketch (1-loop)

Once 18b closes:
- `loops=1` argument to `solve_tree_pipeline` → `solve_loop_pipeline`
- Per-emission propagator denoms → PaVe scalar functions via Layer 5
- Tensor reduction (TID/OPP) for non-trivial loop integrals
- Cross-validation against existing `vertex_g2`, `self_energy_1loop`,
  `running_alpha`, `nlo_box` paths

### What's still deferred to Phase 18b

- **Internal fermion propagators** (Compton tree s+u): walk_fermion_lines
  errors with a deferral message; propagator_num for fermion + composite
  momentum errors deliberately. ~150 LOC to lift.
- **Multi-orbit interference** (Bhabha s+t, multi-channel φ³):
  solve_tree_pipeline currently picks `bundles[1]` as a Phase-18a
  shortcut. Phase 18b sums Burnside-weighted across orbits. ~50 LOC.
- **Boson polarisation** (QCD qq̄→gg with external gluons):
  build_externals returns `(nothing, nothing)` for boson legs.
  ~80 LOC including the polarisation_sum hookup.
- **4-vertex** (gggg): build_vertices errors. ~60 LOC.
- **Symbolic mass** support: external propagator denominators currently
  use `1//1` placeholder for non-zero mass (matches amplitude.jl
  convention). Symbolic mass arrives in 18b. ~40 LOC.
- **Coupling assignment**: AmplitudeBundle.coupling = alg(1) placeholder.

### NEXT SESSION DECISION POINT

Phase 18a closure unlocks several directions:

**Option A — Phase 18b: lift the deferrals (HIGHEST leverage)**
- A1. Multi-orbit Burnside summation in solve_tree_pipeline (~50 LOC)
- A2. Internal fermion propagators (Compton tree validation) (~150 LOC)
- A3. Boson polarisation (QCD qq̄→gg) (~80 LOC)
- A4. Symbolic mass support (~40 LOC)
Estimated 1-2 sessions for a complete tree-level pipeline.

**Option B — 1-loop bridge (Phase 18c)**
After 18b: extend the bridge to 1-loop emissions (PaVe integrals via
existing Layer 5). Phase 18c is the natural sequel and unlocks Spiral 10.

**Option C — Filter ports + golden master push**
nosigma (~120 LOC, +4 cases), floop (~30 LOC, +3 cases) per
Session 24's NEXT SESSION block. Quick wins.

**Recommendation**: A. The 18a milestone proves the bridge works;
18b lifts the artificial scope restrictions to make it useful for
the full tree-level Standard Model.

## SESSION 24 TIMELINE

1. Read `HANDOFF.md`, `Feynfeld_PRD.md`, `CLAUDE.md`. Internalized rules and open bugs.
2. Investigated BUG 2 via GRIND METHOD (read Nogueira 1993, `ALGORITHM.md`,
   qgraf `qgraf-4.0.6.f08:12001-14575`, Julia `src/v2/qgraf/*`).
3. Built `grind/ctrl_phi3_2L.dat` + parsers, ran instrumented qgraf →
   `grind_phi3_2L.txt` (465 emissions, 50 canonical topologies).
4. Ran Julia per-topology dump → 52 topologies, 483 Burnside. Cross-tabbed
   via `compare_topos.jl` → identified 2 excess iso-class forms (+12, +6 = +18).
5. Identified root cause: Julia `step_c_enumerate!` lacks qgraf's post-fill
   permutation canonicality check (f08:13156-13291).
6. Spawned research agent (verified mechanism), implemented `is_canonical_qgraf!`
   in `canonical.jl` + wired into `step_c_enumerate!`. Spawned review agent
   (verified fix, 7/7 criteria).
7. All tests green. Committed as `5e6ddac`. BUG 2 closed.
8. Ran golden master against `count_diagrams_qg21` (the qg21 path):
   loops≤4 → **95 PASS / 0 FAIL / 0 ERROR** (1 initial FAIL on `nosnail` isolated).
9. Traced nosnail discrepancy (25 vs 22) to qgraf `f08:2794-2798`: `nosn>0`
   sets both `intf(nsl)=1` AND `intf(nsb)=1` — nosnail = no-self-loop + no-sbridge.
   Fixed in audition.jl.
10. Phase 17c pipeline swap: `count_diagrams` → `QgrafPort.count_diagrams_qg21`;
    wired all 9 filter kwargs (onepi, nosbridge, notadpole, onshell, nosnail,
    onevi, noselfloop, nodiloop, noparallel).
11. `test_diagram_gen.jl::"QED 1-gen 1-loop"` used `qed_model()` (2-gen) but
    expected qed1's value of 6 for γγ→γγ 1L. Legacy bug was masking this; qg21
    correctly returns 12 for qed2. Fixed to use `qed1_model()`.
12. Full regression: v2 25/25 ✓, qgraf-port 21/21 ✓, golden master 95/104.
    Committed as `d1fa8ee`.

---

## SESSION 24 ACCOMPLISHMENTS — BUG 2 FIXED

**Root cause** (verified by GRIND METHOD with instrumented qgraf-grind in `grind/`):
Julia's `step_c_enumerate!` lacked qgraf's post-fill permutation canonicality
check (qgraf-4.0.6.f08:13156-13291, labels 77/93/102/202/204/114/63). Step C's
cross-row/col checks (f08:12911-12946) are necessary but not sufficient —
without the post-fill perm iteration, Julia emitted 52 canonical topologies
for phi3 2L φφ→φφ (6-deg-3 internal partition) where qgraf emits 50, giving
+18 over-count in `count_diagrams_qg21` (483 vs 465).

The 2 extras were iso-pairs where the rejection perm is INTERNAL-ONLY (not
requiring external swap). qgraf iterates class-respecting perms with `xp(n_ext)`
pinned and rejects when `gam(xp(i1), xp(i2)) > gam(i1, i2)` at any internal
pair (i1 ≥ rhop1). Diagnostic via `grind/compare_topos.jl`:

| Iso-class | qgraf keeps | Julia also kept (extra) | Burnside contrib |
|-----------|------------|------------------------|------------------|
| 1 | gam=(5,6)(5,9)(6,10)(7,8)(7,10)(8,10)(9,9)=2 | gam with (5,10)(6,10)(7,9) | +12 |
| 2 | gam with (7,7)=2(7,8)(9,10)=2 | gam with (7,9)(8,10)=2 | +6 |

Total excess: 12+6 = **18** ✓.

**The fix** (~70 LOC):
- `src/v2/qgraf/canonical.jl`: added `_compare_internal_adjacency(state, perm, rhop1)`
  (internal-only pair comparison, matches f08:13206-13212) and `is_canonical_qgraf!(state)`
  (qgraf-convention lex-LARGEST canonicality, class product with last-ext pinned,
  matches f08:13156-13291).
- `src/v2/qgraf/topology.jl`: `step_c_enumerate!` emit path now calls
  `is_canonical_qgraf!(state)` after `_is_connected_internal` check; reject →
  `@goto row_decrement` for backtrack.
- `src/v2/qgraf/QgrafPort.jl`: export `is_canonical_qgraf!`.

**Verification** (20/20 spot-checks + full test suite):

| Test surface | Result |
|---|---|
| phi3 φφ→φφ 2L (THE bug case) | 483 → **465** ✓ |
| phi3 φ→φφ 2L | 58 ✓ |
| phi3 φφ→φφ 1L | 39 ✓ |
| phi3 φφ→φφ tree | 3 ✓ |
| QED1 15 cases vs golden masters | 15/15 ✓ |
| QED2 ee→μμ 1L | 18 ✓ |
| QCD tree (qq̄→gg, gg→gg, qg→qg) | 3/3 ✓ |
| Full v2 suite (25 files) | 25/25 pass, 0 fail/error |
| qgraf-port tests (21 files) | 21/21 pass (only phase17 B/C still broken) |

**Test-marker updates**:
- `test_qg21_battery.jl`: φ³ 2L 465 case `@test_broken` → `@test`.

**Phase 17c — PIPELINE SWAP COMPLETE**:
- `src/v2/diagram_gen.jl::count_diagrams` now delegates to
  `QgrafPort.count_diagrams_qg21`. Legacy implementation preserved as
  `_count_diagrams_legacy` for regression testing.
- Filter kwargs wired through: `onepi`, `nosbridge`, `notadpole`, `onshell`,
  `nosnail` (= no self-loop + no sbridge per qgraf f08:2794-2798),
  `onevi`, `noselfloop`, `nodiloop`, `noparallel`.
- Golden master coverage jumped from **70/104 → 95/104 PASS** (0 FAIL, 0
  ERROR, 9 SKIP). Remaining SKIPs: 2 qgraf FAIL cases, 4 nosigma, 3 floop.
- Test fix: `test_diagram_gen.jl::"QED 1-gen 1-loop"` now correctly uses
  `qed1_model()` instead of `qed_model()` (2-gen). The previous legacy
  count_diagrams returned 6 for qed2 γγ→γγ 1L (a hidden bug) — the qg21
  path correctly returns 12 (μ-loop included).

**Remaining work on the qg21 port**:
1. Port `nosigma` filter (qgsig, f08:13669) — rejects self-energy insertions.
2. Port `floop` flag (require ≥1 fermion loop).
3. Port `onshellx` (qumvi(3)) and `cycli` filters.
4. Phase 18: Diagram → AlgSum amplitude bridge (the actual payoff — emissions
   carry (xg, ps1, pmap, fermion_sign), convert into Layer 4 AlgSum).

**GRIND diagnostic artefacts** (per-case traces gitignored via
`grind/grind_*.txt`, `grind/julia_*.txt`):
- `grind/ctrl_phi3_2L.dat` — qgraf config for phi3 2L φφ→φφ.
- `grind/parse_grind_trace.jl` — extract per-topology buckets from qgraf trace.
- `grind/dump_julia_phi3_2L.jl` — Julia per-topology Burnside dump.
- `grind/compare_topos.jl` — iso-class cross-tab between qgraf and Julia
  with WL-like signature + per-topology excess detection.

---

## SESSION 23 ACCOMPLISHMENTS — BUG 1 FIXED (preserved for context)

**Root cause** (verified by GRIND METHOD with instrumented qgraf-grind in `grind/`):
qgraf's `dpntro` (rule lookup table built by `qrvi:22020-22090`) stores ALL
distinct positional permutations of each vertex (12 rules for QED2 deg-3, 6 perms
× 2 vertex types). Julia's previous `_qgen_recurse` stored 1 sorted multiset per
fieldset and assigned `_multiset_diff(rule, assigned)` (sorted) to slots in fixed
order, missing emissions where the slot ordering of "remaining" fields was
non-canonical.

For ee→μμ 1L: missing 1 orbit on penguin topology + 1 on box topology
(Burnside contribution 1+1=2; A=16 vs qgraf=18).

**The fix** (`src/v2/qgraf/qgen.jl`, +96/-38 LOC):
- New helper `_qgen_check_perm` implements qgraf's two slot-ordering filters:
  - Self-loop pair check (`qgen:13921-13934`): conjugate pairs in canonical order
  - Multi-edge ordering (`qgen:13948-13954`): consecutive slots to same neighbour sorted
- `_qgen_enumerate_recurse` and `_qgen_recurse` now iterate
  `multiset_permutations(remaining, length(remaining))` per matching multiset rule,
  apply the filters, then recurse if valid.

**Verification** (35/35 spot-checks against qgraf golden masters + full test suite):

| Test surface | Result |
|---|---|
| ee→μμ 1L (THE bug case) | 16 → **18** ✓ |
| QED1 (15 cases vs golden masters) | 15/15 ✓ |
| QED2 (7 cases vs golden masters) | 7/7 ✓ |
| QCD (13 cases vs golden masters) | 13/13 ✓ |
| Phase 17b battery (`test_qg21_battery.jl`) | 23 pass + 1 broken (BUG 2 unchanged) |
| Phase 17 audition (`test_phase17_audition.jl`) | 17 pass + 4 broken (B/C still over-count) |
| Full v2 suite (25 files) | 25/25 pass, 0 fail/error |

**Test-marker updates**:
- `test_count_diagrams_qg21.jl`: ee→μμ 1L `@test_broken` → `@test`
- `test_phase17_audition.jl`: A Burnside `@test_broken` → `@test`; B/C remain
  `@test_broken` (now over-count to 19; the canonicality bug from the
  audition VERDICT is unaffected by Phase 12d).

**Side fix**: `Combinatorics` added to `Project.toml [deps]` (was only transitive
via Manifest; would silently break the build if a dep upgrade dropped the
transitive pull). Reviewer S2.

**BUG 2 status (UNCHANGED — separate root cause)**: φ³ φφ→φφ 2L still returns
483 vs target 465. Phase 12d is provably a no-op for φ³ (single rule, single
multiset perm). The C3 TODO note in `qgen.jl` flags one suspect: `pmap[vv,
rdeg+1..vdeg]` is not saved on backtrack (only neighbour slots are). Benign for
currently-passing cases, but worth checking against BUG 2 where self-loop
topologies abound.

## GRIND METHOD reusable infrastructure

`grind/` directory (qgraf binary/source/traces gitignored):
- `run_v2_tests.sh` — sequential v2 test runner with incremental output
- `dump_julia_emissions.jl` — Julia-side per-emission state dump
- `inspect_dpntro.gdb`, `inspect_qgen.gdb` — gdb scripts for instrumented qgraf
- `ctrl.dat` — qgraf control file for the bug case
- `README.md` — how to instrument and rebuild qgraf locally

Use this same workflow on BUG 2: instrument qgraf, dump per-emission state for
phi3 φφ→φφ 2L (483 vs 465 = +18 over-count), compare against Julia trace,
identify first divergence.

---

## SESSION 22 ACCOMPLISHMENTS

24 commits this session.  The full algorithmic core of qgraf is now ported
to `src/v2/qgraf/` (~1900 LOC, 8 files) with line-by-line citations to
`refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08`.  ~440 tests added.
A `count_diagrams_qg21` entry-point exists as the Strategy C wrapper and
matches legacy on 23 of 25 battery cases.

### Phase-by-phase

| Phase | What | Source ref | Tests | Commit |
|-------|------|------------|------:|--------|
| 2 | Variable-arity `VertexRule` (`NTuple{3}` → `Tuple{Vararg{Symbol}}`) | qgraf model files | 10 | `b0472bb` |
| 3 | 4-gluon vertex `[g,g,g,g]` in `qcd_model` | models/qcd | 3 | `50e90c0` |
| 4 | Faddeev-Popov ghost field + ghost-gluon vertex in `qcd_model` | models/qcd | 3 | `18af6ef` |
| 5 | qg21 Step A audit — degree-seq init in TopoState | f08:12479-12492 | 12 | `fec5b46` |
| 6 | qg21 Step B — `step_b_enumerate!` (xc/xn enumeration) | f08:12554-12658 | 45 | `26763ba` |
| 7a | Step C trivial — single-internal short-circuit | f08:12742-12815 | 6 | `15b1bdb` |
| 7b/c | Step C full state machine (dsum + xg-diag bt + row fill) | f08:12659-13150 | 10 | `16a4f7d` |
| 8 | `_is_connected_internal` BFS at emit point | f08:12980-13038 | 5 | `55196f0` |
| 9 | Integration regression pinning + abstraction invariant | — | 11 | `2bc7552` |
| 11 | qg10 — `qg10_enumerate!` (Knuth Algorithm L) | f08:12001-12200 | 15 | `6cb8b64` |
| 12a | `build_dpntro` lookup table | f08:13889 | 24 | `b6d6fa3` |
| 12b | `compute_qg10_labels` (vlis/vmap/lmap from xg) | f08:12028-12102 | 34 | `2f28165` |
| 12c | `qgen_count_assignments` recursive backtracker | f08:13880-13987 | 3 | `8d23ef5` |
| 13 | `qdis_fermion_sign` (signed half-edge encoding + pair cancellation) | f08:14465-14575 | 1 | `e7aeb3b` |
| 14a | Inline filters: `has_no_selfloop/diloop/parallel` | f08:13960-13978, 13065-13076 | 10 | `5bd43b8` |
| 14b | qumpi family: `is_one_pi`, `has_no_sbridge/tadpole/onshell` | f08:3690-3776 | 4 | `d17df24` |
| 14c | qumvi family: `has_no_snail`, `is_one_vi` | f08:3777-3881 | 6 | `6cfa571` |
| 15 | `compute_local_sym_factor` (S_local) | f08:14361-14411 | 4 | `cc0b00c` |
| 16 | `build_spanning_tree`, `count_chords` | f08:13315-13402 | 5 | `1042da0` |
| 17 prep | `enumerate_topology_automorphisms` (full auto group) | f08:13180-13290 | 13 | `0cff64b` |
| 17a | Dedup audition: 3 strategies + verdict | — | 16+5 broken | `f0ff403`/`bff47ef`/`7842015` |
| 17b | `count_diagrams_qg21` — Strategy C entry point + battery | — | 32+2 broken | `b19d33b` |

### Net golden-master impact

Before session (HANDOFF Session 21, loops ≤ 2):
- PASS: 63 / 104, FAIL: 14, SKIP: 26, ERROR: 1

After session (loops ≤ 2):
- **PASS: 70 / 104** ⬆ +7
- **FAIL: 8** ⬇ −6
- SKIP: 26 (unchanged — Phase 17c needed to wire Phase 14 filters into the legacy pipeline)
- **ERROR: 0** ⬇ −1

Cases moved FAIL/ERROR → PASS (all from phases 2-4 — no pipeline swap needed):
- `qcd gg → gg 0L` (3 → 4)
- `qcd qq̄ → ggg 0L` (15 → 16)
- `qcd gg → ggg 0L` (15 → 25)
- `qcd ghost → ghost 1L onepi` (ERROR → 1)
- `qcd qq̄ → gg 1L onepi` (6 → 7) — gggg-vertex propagation
- `qcd qg → qg 1L onepi` (6 → 7)
- `qcd qq̄ → ggg 0L` (variant)

### Phase 17a audition VERDICT

Three dedup strategies tested against legacy `count_diagrams` on a 10-case
battery.  All three operate on the same emission stream from
`qgen_enumerate_assignments`; disagreement points to the dedup logic.

| Strategy | Score | Diagnosis |
|---|---|---|
| **(A) Burnside** | 9/10 | ✓ **CHOSEN.** `Σ |Stab(emission)| / |G|` over the joint (ps1, pmap) orbit. Robust to in↔out crossings in the auto group. |
| (B) Canonical-pmap | 7/10 | ✗ Compares (ps1, pmap_sig) lex; the canonical orbit-rep may be INVALID for qgen, so the orbit yields 0 emissions instead of 1 (under-count). |
| (C) Pre-filter | 7/10 | ✗ Same bug as (B): pre-filtering ps1 to orbit-reps assumes the rep is qgen-valid, which it may not be. |

**Recommendation**: use (A) Burnside for the Phase 17c pipeline swap.

### Phase 17b: `count_diagrams_qg21` Strategy C entry point

```julia
count_diagrams_qg21(model, in_fields, out_fields; loops=0, onepi=false) -> Int
```

Pipeline: `qg21_enumerate!` → qg10 ext-perm loop → `qgen_enumerate_assignments`
→ Burnside dedup (`Σ |Stab|/|G|`).  Optional `onepi` filter via `is_one_pi`.

**Battery results** (`test/v2/qgraf/test_qg21_battery.jl`): **23 of 25 cases match legacy + qgraf golden master** as integers.  Two outliers documented as `@test_broken` and discussed below.

---

## KNOWN BUGS — BOTH BLOCKERS FIXED

BUG 1 fixed Session 23. BUG 2 fixed Session 24. No blockers remaining for
Phase 17c pipeline swap (only the Phase 14 filter wiring + golden-master
re-verification, per "Pipeline swap (Phase 17c) status" above).

### BUG 1 — qgen flavor-loop under-count (QED multi-gen 1L) — **FIXED Session 23**

**Was**:
```
count_diagrams_qg21(qed_model(), [:e,:e], [:mu,:mu]; loops=1)  →  16  (legacy: 18, qgraf golden: 18)
```

**Now**:
```
count_diagrams_qg21(qed_model(), [:e,:e], [:mu,:mu]; loops=1)  →  18 ✓
```

**Root cause** (verified by GRIND METHOD): qgraf's qgen iterates ALL distinct
positional permutations of each vertex (qrvi:22020-22090); Julia's
`_qgen_recurse` was using multiset matching with sorted slot assignment,
missing valid emissions where the slot ordering of "remaining" fields was
non-canonical. Missed 1 orbit on penguin + 1 on box.

**Fix**: `src/v2/qgraf/qgen.jl` — Phase 12d. New `_qgen_check_perm` helper
implements qgraf's self-loop pair check (qgen:13921-13934) and multi-edge
ordering filter (qgen:13948-13954); `_qgen_{recurse,enumerate_recurse}` now
iterate `multiset_permutations(remaining)` per matching rule and apply the
filters. See SESSION 23 ACCOMPLISHMENTS above.

**Verification**: 35/35 spot-checks against qgraf golden masters
(QED1: 15, QED2: 7, QCD: 13) all match. Full v2 suite: 25/25 pass.

**Diagnostic infra**: see `grind/` (instrumented qgraf gitignored, our
scripts and README committed).

### BUG 2 — phi3 2-loop φφ→φφ over-count (+18) — **FIXED Session 24**

**Was**:
```
count_diagrams_qg21(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=2)  →  483  (legacy: 465, qgraf golden: 465)
```

**Now**:
```
count_diagrams_qg21(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=2)  →  465 ✓
```

**Root cause** (verified by GRIND METHOD): Julia's `step_c_enumerate!` lacked
qgraf's post-fill permutation canonicality check (qgraf-4.0.6.f08:13156-13291).
Step C's cross-row/col checks (f08:12911-12946) are necessary but not sufficient.
Julia emitted 52 canonical topologies for the 6-deg-3 partition where qgraf
emits 50; the 2 extras contributed +12 and +6 to Burnside = **18**.

**Fix**: `src/v2/qgraf/canonical.jl` + `topology.jl`. New function
`is_canonical_qgraf!` implements qgraf's post-fill check: iterates xp over
class product space (externals 1..n_ext-1 + internal (vdeg,xn,xg_diag) classes,
with last-ext pinned), compares internal-pair `gam(xp)` vs `gam(orig)`, rejects
when xp gives lex-LARGER at first-difference position. Called from
`step_c_enumerate!` emit path. See SESSION 24 ACCOMPLISHMENTS above.

**Verification**: 20/20 spot-checks against qgraf golden masters. Full v2
suite 25/25 pass. qgraf-port suite 21/21 pass.

**Diagnostic infra**: see `grind/` — `ctrl_phi3_2L.dat`, `parse_grind_trace.jl`,
`dump_julia_phi3_2L.jl`, `compare_topos.jl`.

### BUG 3 (low priority) — vacuum n_ext=0 not supported by qg10 labels

**Symptom**: `compute_qg10_labels` errors with `qg10_1` ("no candidate vertex with positive vaux") on n_ext=0 partitions (vacuum diagrams).

**Diagnosis**: this matches qgraf's own behavior (qg10:12055-12059 errors identically).  qgraf doesn't generate vacuum diagrams; we don't either.  Vacuum tests are excluded from Phase 12b coverage with a citing comment.

**Where to look**: `src/v2/qgraf/qgen.jl::compute_qg10_labels` — the `vaux=0` branch.  Fix would be to allow lowest-index unvisited at first pick when n_ext=0.

**Priority**: LOW.  Not used by Feynfeld.

---

## CURRENT STATE

### Tests green (in-tree)

| Suite | Count | Status |
|-------|------:|--------|
| `test/v2/test_diagram_gen.jl` | 32 | 32/32 ✓ |
| `test/v2/test_qcd_4gluon.jl` | 3 | 3/3 ✓ (Phase 3) |
| `test/v2/test_qcd_ghost.jl` | 3 | 3/3 ✓ (Phase 4) |
| `test/v2/test_vertex_arity.jl` | 10 | 10/10 ✓ (Phase 2) |
| `test/v2/qgraf/test_types.jl` | 98 | 98/98 ✓ |
| `test/v2/qgraf/test_canonical.jl` | 35 | 35/35 ✓ |
| `test/v2/qgraf/test_step_a.jl` | 12 | 12/12 ✓ (Phase 5) |
| `test/v2/qgraf/test_step_b.jl` | 45 | 45/45 ✓ (Phase 6) |
| `test/v2/qgraf/test_step_c.jl` | 10 | 10/10 ✓ (Phase 7) |
| `test/v2/qgraf/test_step_c_connectedness.jl` | 5 | 5/5 ✓ (Phase 8) |
| `test/v2/qgraf/test_qg21_integration.jl` | 11 | 11/11 ✓ (Phase 9) |
| `test/v2/qgraf/test_qg10.jl` | 15 | 15/15 ✓ (Phase 11) |
| `test/v2/qgraf/test_qgen_dpntro.jl` | 24 | 24/24 ✓ (Phase 12a) |
| `test/v2/qgraf/test_qg10_labels.jl` | 34 | 34/34 ✓ (Phase 12b) |
| `test/v2/qgraf/test_qgen_recurse.jl` | 3 | 3/3 ✓ (Phase 12c) |
| `test/v2/qgraf/test_qdis.jl` | 1 | 1/1 ✓ (Phase 13) |
| `test/v2/qgraf/test_filters_inline.jl` | 10 | 10/10 ✓ (Phase 14a) |
| `test/v2/qgraf/test_filters_qumpi.jl` | 4 | 4/4 ✓ (Phase 14b) |
| `test/v2/qgraf/test_filters_qumvi.jl` | 6 | 6/6 ✓ (Phase 14c) |
| `test/v2/qgraf/test_sym_factor.jl` | 4 | 4/4 ✓ (Phase 15) |
| `test/v2/qgraf/test_momentum.jl` | 5 | 5/5 ✓ (Phase 16) |
| `test/v2/qgraf/test_automorphisms.jl` | 13 | 13/13 ✓ (Phase 17 prep) |
| `test/v2/qgraf/test_phase17_audition.jl` | 17+4 broken | (Phase 17a — B/C dedup still broken) |
| `test/v2/qgraf/test_count_diagrams_qg21.jl` | 10 | 10/10 ✓ (Phase 17b) |
| `test/v2/qgraf/test_qg21_battery.jl` | 24 | 24/24 ✓ (Phase 17b, BUG 2 fixed Session 24) |

Plus full v2 regression (test_ee_mumu_x, test_ee_ww, test_qqbar_gg, etc.) all green.

### Files in `src/v2/qgraf/`

| Path | LOC | Purpose |
|------|----:|---------|
| `QgrafPort.jl` | 25 | submodule wrapper + exports |
| `types.jl` | 195 | Partition, EquivClass, FilterSet, TopoState |
| `canonical.jl` | ~360 | is_canonical_full!, is_canonical_qgraf! (Session 24), enumerate_topology_automorphisms |
| `topology.jl` | ~520 | step_b_enumerate!, step_c_enumerate!, qg10_enumerate!, _is_connected_internal |
| `qgen.jl` | ~330 | build_dpntro, compute_qg10_labels, qgen_count_assignments, qgen_enumerate_assignments, qdis_fermion_sign, compute_local_sym_factor |
| `filters.jl` | ~190 | has_no_*, is_one_pi, is_one_vi |
| `momentum.jl` | ~80 | build_spanning_tree, count_chords |
| `audition.jl` | ~290 | count_dedup_burnside/canonical/prefilter, count_diagrams_qg21, is_emission_canonical, emission_stabilizer |

### Files outside `src/v2/qgraf/` modified

| Path | What |
|------|------|
| `src/v2/rules.jl` | VertexRule.fields → `Tuple{Vararg{Symbol}}`; FeynmanRules.vertices → `Dict{Tuple, VertexRule}` |
| `src/v2/qcd_model.jl` | gggg vertex + ghost field |
| `src/v2/ew_model.jl` | Dict signature relaxed |
| `src/v2/phi3_model.jl` | Dict signature relaxed |
| `test/v2/test_diagram_gen.jl` | gg→gg test asserts 4 (was 3, documented as known gap) |

---

## REMAINING WORK

### Phase 17c — pipeline swap (gated on BUG 2)

Replace `count_diagrams` in `src/v2/diagram_gen.jl` with a wrapper that
calls `count_diagrams_qg21`.  Wire the Phase 14 filter predicates into the
new path so the 26 SKIP cases unlock.

**Gating**: BUG 1 fixed Session 23. BUG 2 (φ³ 2L φφ→φφ over-count +18)
still gates the swap; otherwise Phase 17c regresses test_diagram_gen on
the φ³ 2L 4-point case.

**BUG 2 next-step diagnostic** (after fixing it): re-run
`grind/run_v2_tests.sh` and the 35 spot-check battery; rerun the golden
master report to see how many of the 26 currently-SKIP cases turn green.

### Phase 18 — Diagram → AlgSum amplitude bridge (Layer 4)

Each emission from the new pipeline carries (xg, ps1, pmap, fermion_sign).
Convert this into the existing AlgSum amplitude structure used by the
v2 algebra layer.  Required for actual amplitude evaluation.

### Other deferred work

- qpg11 partition iterator (Phase 10): currently uses legacy
  `_degree_partitions`; works fine.  Could port faithfully later.
- Filter integration into `count_diagrams_qg21`: current API only handles
  `onepi`; extend to the full FilterSet.
- `qgsig` (nosigma) and `qcyc` (cycli): need momentum routing first
  (they consume qgraf's `flow[][]` array — Phase 16 deferred work).
- Full S_nonlocal: extend `enumerate_topology_automorphisms` to include
  ext-perm orbits that preserve the field assignment (currently we use
  topology-only autos and let Burnside handle it).

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. Critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **CITE EVERYTHING.** Local file path + equation number + verbatim equation.
3. **ALL TESTS SYMBOLIC.** No numerical spot-checks. AlgSum == AlgSum only.
4. **JULIA IDIOMATIC.** Dispatch, not isa cascades. No Any.
5. **NO PARALLEL JULIA AGENTS.** Read-only research CAN run in parallel.
6. **LOC LIMIT ~200.** No source file exceeds ~200 lines.
7. **REVIEW.** Rigorous reviewer after every core change.
8. **TIERED WORKFLOW.** Core (>20 LOC): 3 research + 1 review.

---

## QUICK COMMANDS

```bash
# Full qgraf-port test suite
for f in test/v2/qgraf/*.jl; do julia --project=. "$f"; done

# Main regression
julia --project=. test/v2/test_diagram_gen.jl            # 32/32

# New Strategy C entry point — quick sanity
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
using .FeynfeldX.QgrafPort: count_diagrams_qg21
println(count_diagrams_qg21(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=1))   # 39
println(count_diagrams_qg21(qed_model(), [:e,:e], [:mu,:mu]; loops=0))           # 1
'

# Audition battery vs legacy (10 cases, A/B/C comparison)
julia --project=. scripts/audition_compare.jl

# Golden-master diagnostic
QGRAF_MAX_SECONDS=60 julia --project=. scripts/qgraf_golden_master_report.jl 2

# Beads
bd ready
bd stats
bd dolt push

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```

---

## SESSION 21 CONTEXT (preserved)

The 474 → 465 phi3 2-loop canonicality fix from Session 21 is INTACT.
`topology_filter._is_canonical_topo` still delegates to
`QgrafPort.is_canonical_feynman`.  Verified by `test/v2/qgraf/test_qg21_integration.jl`
asserting `legacy_count(phi3, [:phi,:phi], [:phi,:phi], loops=2) == 465`.

Note: BUG 2 (fixed Session 24) was on the `count_diagrams_qg21` path (returned
483), not the legacy path (always returned 465 correctly). Post Phase 17c,
`count_diagrams` now delegates to `count_diagrams_qg21`, so both paths agree
and reach 465 for this case.

---

## NEXT SESSION DECISION POINT

Both BUG 1 and BUG 2 are fixed. Phase 17c pipeline swap is complete. The
qg21 port is the default counting path. Three candidate directions for
the next session, in rough order of payoff vs effort:

### Option A — Phase 18: Diagram → AlgSum amplitude bridge (HIGHEST leverage)

This is where the qg21 port actually **pays off**. Currently the pipeline
only COUNTS diagrams; it doesn't produce amplitudes. Phase 18 bridges the
gap: each emission `(xg, ps1, pmap, fermion_sign)` becomes an `AlgSum` in
Layer 4, which Layer 5 (PaVe reduction) and Layer 6 (cross section) consume.

**Estimated effort**: ~600 LOC total, ~2-3 sessions.

| Subtask | LOC | Notes |
|---|---|---|
| A1. Complete momentum routing (Phase 16 partial) | ~120 | Spanning tree + leaf-peeling exist; need per-edge momentum assignment + sign normalization. Citation: ALGORITHM.md §5.1-5.2, qgraf f08 `flow[][]` array. |
| A2. Emission → AlgSum builder | ~250 | For each emission, construct: propagator factors × vertex factors × ext spinors/pol × fermion sign × 1/S prefactor. |
| A3. Wire into `cross_section.jl` | ~80 | Replace hand-built `tree_channels` / `loop_channels` with pipeline-generated input. |
| A4. Validation tests | ~150 | Reproduce ee→μμ tree + 1L, Compton, Bhabha via pipeline, cross-check against existing hand-built AlgSums to machine precision. |

**Risk factors (could 2x the estimate)**:
1. Layer 4 AlgSum may need small extensions to accept pipeline-shaped inputs.
2. `qdis_fermion_sign` returns only ±1; Phase 18 needs the full trace
   ordering (directed traversal of fermion lines).
3. Symmetry factor `1/S` — `S_local` exists (`qgen.jl::compute_local_sym_factor`),
   but `S_nonlocal` is currently computed aggregate-only via Burnside
   `|Stab|/|G|`; Phase 18 needs it per-emission.
4. 1-loop cases force Layer 5 (PaVe) interaction — can defer as Phase 18b.

**Suggested scoping**: **Phase 18a = tree-level only** first (~300 LOC,
1 session): momentum routing + AlgSum builder + ee→μμ-tree validation.
Defer 1-loop to Phase 18b. Visible physics payoff, bounded risk.

### Option B — port remaining filters (modest payoff, small LOC)

The 9 golden-master SKIPs break down as:
- **2 qgraf FAIL cases** (not fixable — qgraf itself can't generate them).
- **4 `nosigma` cases** — `qgsig` at qgraf f08:13669 rejects self-energy
  insertions. Requires BFS-based 2-point subdiagram detection. **~80-120 LOC.**
- **3 `floop` cases** — fermion-loop counter + filter. Infrastructure
  partially exists in `qgen.jl` (`antiq` tracking at f08:13988-14034).
  Expose count + compare. **~30 LOC.**

`floop` is cheap and unlocks 3 cases. `nosigma` is moderate and unlocks 4.
Pure counter-mode improvements; no new physics capability.

### Option C — Spiral 8 remainder (chiral physics unblock)

- γ5 traces (`feynfeld-qu1`): unblocks chiral EW.
- Eps (Levi-Civita) contraction completion.
- MUnit translation continues alongside (per revised PRD §3.3).

This is Layer 4 work, independent of the qg21 port. Parallel track —
could be picked up by any agent that has capacity.

### Recommendation

**Option A (Phase 18a — tree-level)**: highest payoff. The qg21 port is
a diagram counter that doesn't do physics yet. Phase 18a ends with the
pipeline producing ee→μμ tree-level amplitudes matching the existing
hand-built implementation — a visible, bounded milestone that the rest
of the architecture (Layers 5, 6) can consume.

If the next session prefers a quick win first, knock off `floop` (~30
LOC, unlocks 3 golden masters) as a warm-up, then Phase 18a.
