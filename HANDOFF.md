# HANDOFF — 2026-04-14 (Session 24: BUG 2 + Phase 17c pipeline swap)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms.
2. Run `bd ready` to see available work.
3. Run `julia --project=. test/v2/test_diagram_gen.jl` → expect 32/32 green.
4. Run all qgraf-port tests: `for f in test/v2/qgraf/*.jl; do julia --project=. "$f"; done`
   → all green (only phase17 audition has `@test_broken` markers for B/C
   dedup strategies — known, unrelated to BUG 1/2).
5. Run full v2 suite: `grind/run_v2_tests.sh` → 25/25 files pass, 0 fail/error.
6. Optional: `QGRAF_MAX_SECONDS=120 julia --project=. scripts/qgraf_golden_master_report.jl 2`
   → legacy `count_diagrams` path still PASS=70/104 (unchanged by BUG 2 fix).
   `count_diagrams_qg21` path now matches qgraf on all tested cases.

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

Note: BUG 2 above is on the NEW `count_diagrams_qg21` path (returns 483),
not the legacy path (still returns 465 correctly).  The Session 21 fix is
specific to the legacy enumerator and unaffected by the Strategy C work.
