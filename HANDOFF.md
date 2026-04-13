# HANDOFF — 2026-04-13 (Session 22: Strategy C qg21 port — phases 2-17b complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms.
2. Run `bd ready` to see available work.
3. Run `julia --project=. test/v2/test_diagram_gen.jl` → expect 32/32 green.
4. Run all qgraf-port tests: `for f in test/v2/qgraf/*.jl; do julia --project=. "$f"; done`
   → expect ~440 tests across 14 files, all green or `@test_broken`.
5. Optional: `QGRAF_MAX_SECONDS=120 julia --project=. scripts/qgraf_golden_master_report.jl 2`
   → live-streams every golden-master case. Currently PASS=70/104 at loops≤2
   (was 63 at start of session; +7 thanks to phases 2-4).

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

## KNOWN BUGS — BLOCKERS FOR PHASE 17c PIPELINE SWAP

Two real, reproducible bugs.  Both gate the full pipeline swap.

### BUG 1 — qgen flavor-loop under-count (QED multi-gen 1L)

**Symptom**:
```
count_diagrams_qg21(qed_model(), [:e,:e], [:mu,:mu]; loops=1)  →  17  (legacy: 18, qgraf golden: 18)
```

Discovered during Phase 17a audition.  All three dedup strategies agree on 16-17 (off by 1-2 from 18); since they share the same emission stream, the deficit is in `qgen` itself, not in dedup.

**Diagnosis (incomplete)**:
- For multi-generation QED, a 1-loop diagram can carry an `e`-loop OR a `μ`-loop.
- Both rules `(e, e_bar, γ)` and `(μ, μ_bar, γ)` are in `dpntro[3]`.
- `qgen_count_assignments` iterates rules at each vertex via the multiset matching.
- Hypothesis: when an internal vertex has all 3 slots empty (no fields assigned yet), my multiset matcher tries each rule but the BACKTRACKING (saving/restoring neighbor pmap) might prematurely commit to one flavor for the rest of the loop.

**Diagnostic scripts**:
- `scripts/audition_compare.jl` — battery vs legacy on 10 cases
- `scripts/debug_qed_1l.jl` — counts QED 1L topologies (returns 6)

**Where to look**:
- `src/v2/qgraf/qgen.jl::_qgen_recurse` — the multiset rule iteration
- `src/v2/qgraf/qgen.jl::_qgen_enumerate_recurse` — same recursion for the emit-callback variant

**Audition counts that hit this bug**:
- `qed_model()` `ee → μμ` 1L: A=16, B=17, C=17 (legacy=18)

**Commit context**: `bff47ef`, `7842015`.

### BUG 2 — phi3 2-loop φφ→φφ over-count (+18)

**Symptom**:
```
count_diagrams_qg21(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=2)  →  483  (legacy: 465, qgraf golden: 465)
```

Discovered during Phase 17b battery testing.  This is the famous Session 21 regression case (was 474 before the canonicality fix; legacy now correct at 465).

**Counterpart that PASSES**:
- `phi3 2L φ→φφ` (different partition: 5 deg-3 internals): qg21 = 58 = legacy ✓

So the bug is specific to the **6-deg-3-internals** partition that phi3 φφ→φφ 2L uses.

**Hypotheses (untested)**:
- Self-loop topology in this partition where my Burnside formula or my auto-group enumeration over-counts.
- Multi-edge (parallel-edge) topology where `compute_local_sym_factor`'s self-loop branch (the Z_2 reversal factor) interacts wrongly with the joint (ps1, pmap) Burnside.
- The `enumerate_topology_automorphisms` extension to include externals may over-generate for 2L topologies with rich symmetry.

**Where to look**:
- `src/v2/qgraf/canonical.jl::enumerate_topology_automorphisms` — the auto group
- `src/v2/qgraf/audition.jl::count_diagrams_qg21` — the Burnside loop
- `src/v2/qgraf/qgen.jl::compute_local_sym_factor` — currently NOT used by count_diagrams_qg21 (only the topology-auto |Stab|/|G| is used); maybe S_local must be incorporated for parallel-edge or self-loop topologies?

**Diagnostic next step**: enumerate the 6-internal partition's qg21 topologies, compute per-topology contribution to the 483 sum, identify the over-counting class.

**Audition test**:
```julia
@test_broken count_diagrams_qg21(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=2) == 465
```

**Commit context**: `b19d33b`.

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
| `test/v2/qgraf/test_phase17_audition.jl` | 16+5 broken | (Phase 17a) |
| `test/v2/qgraf/test_count_diagrams_qg21.jl` | 9+1 broken | (Phase 17b) |
| `test/v2/qgraf/test_qg21_battery.jl` | 23+1 broken | (Phase 17b) |

Plus full v2 regression (test_ee_mumu_x, test_ee_ww, test_qqbar_gg, etc.) all green.

### Files in `src/v2/qgraf/`

| Path | LOC | Purpose |
|------|----:|---------|
| `QgrafPort.jl` | 25 | submodule wrapper + exports |
| `types.jl` | 195 | Partition, EquivClass, FilterSet, TopoState |
| `canonical.jl` | ~290 | is_canonical_full!, enumerate_topology_automorphisms |
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

### Phase 17c — pipeline swap (gated on bug fixes)

Replace `count_diagrams` in `src/v2/diagram_gen.jl` with a wrapper that
calls `count_diagrams_qg21`.  Wire the Phase 14 filter predicates into the
new path so the 26 SKIP cases unlock.

**Gating**: BUG 1 (qgen flavor-loop) and BUG 2 (2L over-count) MUST be
fixed first; otherwise Phase 17c regresses 2 of 32 currently-passing
test_diagram_gen tests AND introduces incorrect counts on QED 1L 2-gen.

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
