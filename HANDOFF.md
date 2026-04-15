# HANDOFF — 2026-04-15 (Session 25: Phase 18a — Diagram → AlgSum bridge)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms.
2. Run `bd ready` to see available work.
3. Run `julia --project=. test/v2/test_diagram_gen.jl` → expect all green.
4. Run all qgraf-port tests: `for f in test/v2/qgraf/*.jl; do julia --project=. "$f"; done`
   → 30/30 files green (only `test_phase17_audition.jl` has `@test_broken`
   markers for Strategy-B/C dedup — known audition-verdict issue, unrelated).
5. Run full v2 suite: `./grind/run_v2_tests.sh` → 26/26 files pass, 0 fail/error.
6. Golden master report (loops ≤ 4):
   - `QGRAF_MAX_SECONDS=300 julia --project=. scripts/qgraf_golden_master_report.jl 4`
     → **PASS=95 / FAIL=0 / SKIP=9 / ERROR=0** (unchanged from Session 24).
7. **NEW**: Pipeline ≡ handbuilt validated for ee→μμ tree:
   `julia --project=. test/v2/qgraf/test_phase18a_pipeline.jl` → 1/1 ✓.
8. See **"NEXT SESSION DECISION POINT"** at the bottom for Phase 18b options.

---

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
