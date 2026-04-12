# HANDOFF — 2026-04-12 (Session 21: qgraf Strategy C port — 474→465 fixed)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms.
2. Run `bd ready` to see available work (57 open issues; see §"Open beads" below for the new ones from this session).
3. Run `julia --project=. test/v2/test_diagram_gen.jl` → expect 32/32 green.
4. Run `julia --project=. test/v2/qgraf/test_types.jl` → expect 98/98.
5. Run `julia --project=. test/v2/qgraf/test_canonical.jl` → expect 35/35.
6. Optional: `QGRAF_MAX_SECONDS=120 julia --project=. scripts/qgraf_golden_master_report.jl 2`
   → live-streams every golden-master case with status + timing. Currently PASS=63
   of 104 at loops≤2 (14 FAIL, 26 SKIP, 1 ERROR; breakdown in §"Remaining work").

---

## SESSION 21 ACCOMPLISHMENTS

### 1. **phi3 φφ→φφ 2-loop regression fixed (474 → 465)**

The centerpiece of HANDOFF Session 20's CRITICAL FAILURE section. Root cause was
`topology_enum.jl`'s `_is_canonical_topo` using pairwise swaps over degree-only
equivalence classes — misses multi-vertex automorphisms, producing 9 duplicate
topologies at 2-loop+.

Fix: full per-equivalence-class lex-next-permutation canonicality check ported
from qgraf-4.0.6.f08:13156–13291 (labels 77/93/102/202/204).

`topology_filter._is_canonical_topo` now delegates to
`QgrafPort.is_canonical_feynman` which:
- builds degree-only classes (matching `topology_enum.jl`'s non-sorted fill),
- enumerates full S_n permutations per class (not pairwise swaps),
- compares ALL vertex pairs (ext-int included since no Step A/B pre-sorting),
- uses lex-smallest convention to match the descending fill order.

Commit: `182a657` "Fix qg21 canonicality: 474 → 465 at φ³ 2-loop".

### 2. **Strategy C qgraf port — foundation delivered**

After a 3-way design audition (A: faithful coroutine, B: algebraic pipe,
C: hybrid recursive-descent), picked Strategy C: plain `for` loops outside,
recursive backtracking inside, zero-allocation hot path via preallocated
`TopoState`, callback monomorphisation via `where {F}`.

New submodule `src/v2/qgraf/` under `FeynfeldX.QgrafPort`:

| File | LOC | Contents |
|------|-----|----------|
| `QgrafPort.jl` | 16 | submodule wrapper |
| `types.jl` | 190 | `Partition`, `EquivClass`, `FilterSet`, `TopoState` + helpers; `MAX_V=24` |
| `canonical.jl` | 212 | `compute_equiv_classes!`, `_lex_next_class!`, `next_class_perm!`, `is_canonical_full!`, `is_canonical_feynman` |

Ground truth acquired: **`refs/papers/Nogueira1993_JCompPhys105_279.pdf`**
(fetched via Playwright + TIB VPN, committed to gitignored `refs/papers/`).

### 3. **Multi-generation QED models**

`QEDModel` refactored to hold `Vector{Field{Fermion}}` (was hard-coded electron
+ muon). New constructors:

- `qed1_model()` — electron only, matches qgraf qed1
- `qed_model()` — 2-gen (electron + muon), matches qgraf qed2, unchanged default
- `qed3_model()` — 3-gen (e, μ, τ), matches qgraf qed3

Legacy `.electron`/`.muon`/`.tau` accessors preserved via `getproperty`.

Commit: `349e3a1` "Multi-gen QED models + golden-master diagnostic".

### 4. **Golden-master diagnostic script**

`scripts/qgraf_golden_master_report.jl` parses `SUMMARY.md` (zero deps),
streams per-case PASS/FAIL/SKIP/SLOW/ERROR status live, supports
`QGRAF_MAX_SECONDS` env var for per-case budget and CLI `max_loops` arg.

### 5. **Beads remote configured**

Beads Dolt had no remote configured (bd dolt push was failing silently).
Added `origin → file:///home/tobiasosborne/Projects/Feynfeld.jl/.beads-remote`
as a local backup. `bd dolt push` now succeeds.

Note: beads issue data also persists via `.beads/backup/backup_state.json`
(committed to git) even without a Dolt remote.

---

## CURRENT STATE

### Tests green

| Suite | Count | Status |
|-------|------:|--------|
| `test/v2/test_diagram_gen.jl` | 32 | **32/32 green** (was 31 + 1 RED target) |
| `test/v2/qgraf/test_types.jl` | 98 | 98/98 green (new this session) |
| `test/v2/qgraf/test_canonical.jl` | 35 | 35/35 green (new this session) |
| `test/v2/test_box_ee_mumu.jl` | 136 | 136/136 green (regression check) |
| Other v2 tests | various | spot-checked green (bhabha, self-energy, coeff, pave) |

### Golden-master landscape (loops ≤ 2, 104 cases)

- **PASS: 63** including phi3 φφ→φφ 2L=465, phi3 φ→φφ 2L=58,
  qed1 e→eγ 2L=50, qed1 γ→ee 2L=50, all QED 2-gen onepi, phi3 all tree+1L
- **FAIL: 14** — categorised below
- **SKIP: 26** — cases using filters we haven't ported yet (notadpole/nosnail/nosigma/etc.)
- **ERROR: 1** — qcd ghost→ghost (ghost field not in our qcd_model)

### Files in play

| Path | Purpose | Status |
|------|---------|--------|
| `src/v2/qgraf/QgrafPort.jl` | new submodule | Phase 1 green |
| `src/v2/qgraf/types.jl` | core types | Phase 1 green |
| `src/v2/qgraf/canonical.jl` | per-class lex-perm check | Phase 2 green |
| `src/v2/topology_enum.jl` | legacy enum | delegates canonicality to QgrafPort |
| `src/v2/topology_filter.jl` | legacy filters | `_is_canonical_topo` delegates |
| `src/v2/model.jl` | model types | `QEDModel` refactored |
| `scripts/qgraf_golden_master_report.jl` | diagnostic | NEW |
| `refs/papers/Nogueira1993_JCompPhys105_279.pdf` | canonical algorithm paper | NEW |
| `refs/papers/fetch_nogueira.mjs` | playwright fetch script | NEW |

---

## REMAINING WORK (failures taxonomy)

### A. QCD 4-gluon vertex + ghost field — 6 failing cases

| Case | Got | Want |
|---|---:|---:|
| qcd gg→gg 0L | 3 | 4 |
| qcd qq̄→ggg 0L | 15 | 16 |
| qcd gg→ggg 0L | 15 | 25 |
| qcd gg→gg 1L onepi | 6 | 24 |
| qcd qq̄→gg 1L onepi | 6 | 7 |
| qcd qg→qg 1L onepi | 6 | 7 |
| qcd ghost→ghost 1L onepi | ERROR | 1 |

Structural: `VertexRule.fields::NTuple{3,Symbol}` is hard-coded 3-point.
Needs generalisation to variable arity. Beads: **feynfeld-dr6**.

Cleanest route: change `VertexRule.fields` to `Tuple{Vararg{Symbol}}` and
`FeynmanRules.vertices` to `Dict{Tuple, VertexRule}`. Audit call sites.
`_model_vertex_degrees` and `_degree_partitions` already support arbitrary
degrees. `ExpandedModel.vertex_rules` already uses `Set{Vector{Symbol}}`.

### B. QED fermion-orientation counting — 5-7 failing cases

Current state (after Phase 6 + Phase 2 fixes):
| Case | Got | Want |
|---|---:|---:|
| qed1 γγ→γγ 1L | 3 | 6 |
| qed1 ee→γγ 1L | 23 | 24 |
| qed1 eγ→eγ 1L | 23 | 24 |
| qed1 γγ→ee 1L | 23 | 24 |
| qed1 γ→γ 2L | 4 | 6 |
| qed1 e→e 2L | 11 | 10 |
| qed3 γ→γ 2L onepi | 6 | 9 |

Root cause: `field_assign.jl::_count_closed_loops_expanded` applies
`÷ 2^n_fl` halving per "closed fermion loop". This is partially correct:
- bubbles (multi-edge) already handled by canonical parallel-edge ordering
  (halving skipped via multi-edge check) → ✓
- boxes without external fermion attachment (γγ→γγ) are halved but should not
  be (each fermion-flow orientation is a distinct qgraf diagram) → ✗

Naive removal of halving broke γ→γ 2L, e→e 2L, and ee→γγ 1L (which rely on
the halving for some sub-diagrams). I attempted it in this session and reverted.

**Proper fix**: port qgraf's **qdis** (fermion sign + orientation labelling
via signed half-edge integers) as a separate algorithmic pass, replacing the
heuristic halving in `_count_closed_loops_expanded`. See qgraf-4.0.6.f08
lines 14465–14575. Beads: create new issue.

### C. Filter predicates not implemented — 26 SKIP cases

`count_diagrams` supports only the `onepi` flag. Golden masters use
`notadpole, nosnail, nosigma, noselfloop, nodiloop, noparallel, onevi, floop`.

Each maps to a qgraf filter:
- `qumpi(situ)` dispatch (qgraf-4.0.6.f08:3690) → notadpole/nosbridge/onshell
- `qumvi(situ)` dispatch (qgraf-4.0.6.f08:3777) → nosnail/onevi/onshellx
- `qgsig` (qgraf-4.0.6.f08:13669) → nosigma
- `qcyc` (qgraf-4.0.6.f08:18830) → cycli
- `noselfloop/nodiloop/noparallel` — inline xg-diagonal/entry bounds in Step C

Beads: **feynfeld-con**. `FilterSet` type already defined in
`src/v2/qgraf/types.jl`.

### D. Full Strategy C field-assignment port

Eventually the legacy `field_assign.jl` + `vertex_check.jl` should be
replaced with a Strategy C port of qgraf's `qgen` (dpntro backtracker)
+ `qdis` (fermion sign). This would solve (B) definitively and integrate
with the new `TopoState`/`EquivClass` types.

Design in Session 21 audition document; Strategy C §7 in the design doc
spells it out. Beads: **feynfeld-ney** (master).

---

## OPEN BEADS (Session 21)

- `feynfeld-ney` — P0 master: qgraf Strategy C port (hybrid recursive-descent qg21)
- `feynfeld-chm` — P1: Multi-gen QED models — **DONE, close with --force**
- `feynfeld-dr6` — P1: QCD 4-gluon vertex (gggg) + ghost
- `feynfeld-con` — P1: Filter predicates (notadpole/nosnail/nosigma/etc.)
- `feynfeld-d82` — P1: QED photon `external=` + `notadpole=` keywords
  (hypothesis disproven this session; close with reason)

Closed this session: `feynfeld-nwe` (Phase 0 harness), `feynfeld-5hr`
(Phase 1 skeleton), `feynfeld-xjc` (Phase 2 canonical), `feynfeld-4tg`
(THE 474→465 regression — fixed).

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
# Main regression test
julia --project=. test/v2/test_diagram_gen.jl            # expect 32/32

# New qgraf port unit tests
julia --project=. test/v2/qgraf/test_types.jl            # expect 98/98
julia --project=. test/v2/qgraf/test_canonical.jl        # expect 35/35

# Golden-master diagnostic (LIVE streaming, one line per case)
QGRAF_MAX_SECONDS=120 julia --project=. scripts/qgraf_golden_master_report.jl 2

# Quick count checks
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
println(count_diagrams(phi3_model(), [:phi,:phi], [:phi,:phi]; loops=2))  # 465 (was 474)
println(count_diagrams(qed_model(), [:e,:e], [:mu,:mu]; loops=1))          # 18
'

# Beads
bd ready                                      # available work
bd stats                                       # project health (218 total, 57 open)
bd show feynfeld-ney                          # Strategy C master issue
bd dolt push                                   # now works (local file:// remote)

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```

---

## COMMIT HISTORY (Session 21)

- `182a657` Fix qg21 canonicality: 474 → 465 at φ³ 2-loop (full per-class perm)
- `349e3a1` Multi-gen QED models + golden-master diagnostic
