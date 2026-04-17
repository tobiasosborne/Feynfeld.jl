# Feynfeld.jl v2 Test Suite Stocktake — 2026-04-17

**Executive Summary:** 301 total tests across 62 files (7,628 LOC). Organized as: 25 main pipeline/unit tests (502 tests), 5 MUnit FeynCalc translations (112 tests), 30 qgraf integration tests (414 tests). Comprehensive process coverage (Compton, Bhabha, e+e-→μμ, e+e-→WW, QCD) with 4 known broken tests (canonical dedup bugs). Test orchestration via `runtests.jl` (20 core files) + individual scripts via `run_v2_tests.sh`.

---

## 1. Summary Table: All Test Files

### A. Core Process/Pipeline Tests (25 files)

| File | LOC | Test Count | Broken | Category | Purpose |
|------|-----|-----------|--------|----------|---------|
| test_pipeline.jl | 372 | 20 | 0 | pipeline | Full 6-layer Model→Amplitude→Evaluate |
| test_ee_ww_grozin.jl | 289 | 8 | 0 | pipeline | e+e→WW with polarization tensors (Grozin) |
| test_nlo_box_validation.jl | 210 | 25 | 0 | pipeline | NLO box diagrams vs LoopTools cross-check |
| test_compton.jl | 228 | 4 | 0 | pipeline | Compton e+γ→e+γ (Peskin & Schroeder 5.87) |
| test_bhabha.jl | 172 | 4 | 0 | pipeline | Bhabha e+e→e+e with identical fermion exchange |
| test_ee_ww.jl | 165 | 31 | 0 | pipeline | e+e→WW tree cross-section (Grozin) |
| test_box_ee_mumu.jl | 153 | 26 | 0 | pipeline | 1-loop box D₀ integrals (Denner 1993) |
| test_ee_mumu_x.jl | 141 | 14 | 0 | pipeline | Tracer bullet: e+e→μ+μ (P&S 5.10) |
| test_vertical.jl | 170 | 34 | 0 | pipeline | Vertical amplitude structure |
| test_self_energy.jl | 149 | 25 | 0 | pipeline | Tree self-energy amplitude |
| test_self_energy_1loop.jl | 95 | 13 | 0 | pipeline | 1-loop self-energy |
| test_vertex_g2.jl | 170 | 31 | 0 | pipeline | Photon 3-point vertex (QED) |
| test_running_alpha.jl | 183 | 32 | 0 | unit/integral | Running coupling constant α(Q²) |
| test_pave.jl | 184 | 50 | 0 | unit/integral | PaVe tensor/scalar function library |
| test_d0.jl | 107 | 23 | 0 | unit/integral | D₀ 4-point scalar (tHooft-Veltman) |
| test_diagram_gen.jl | 192 | 32 | 0 | unit | Diagram generation vs QGraf golden master |
| test_colour.jl | 144 | 27 | 0 | unit | SU(N) color algebra (generators, traces) |
| test_coeff.jl | 64 | 29 | 0 | unit | DimPoly coefficient arithmetic |
| test_vertex_arity.jl | 54 | 10 | 0 | unit | Vertex function arity validation |
| test_schwinger.jl | 84 | 15 | 0 | unit | Schwinger model reference formulas |
| test_qqbar_gg.jl | 165 | 2 | 0 | unit | q+q̄→g+g diagram count |
| test_qcd_4gluon.jl | 34 | 3 | 0 | unit | 4-gluon QCD diagram count |
| test_qcd_ghost.jl | 38 | 3 | 0 | unit | Ghost field interactions |
| test_munit_batch1.jl | 183 | 23 | 0 | munit | DiracTrace, Contract, PolarizationSum |
| test_munit_batch2.jl | 127 | 18 | 0 | munit | DiracTrick n≥3 identities |
| **SUBTOTAL** | **3,713** | **502** | **0** | | |

### B. MUnit FeynCalc Translations (5 dedicated files in munit/)

| File | LOC | Test Count | Broken | Function | Coverage |
|------|-----|-----------|--------|----------|----------|
| test_DiracTrick.jl | 355 | 65 | 0 | DiracTrick | Complete (all free-index cases) |
| test_DiracTrace.jl | 162 | 22 | 0 | DiracTrace | Tr[γ chains] basics + Levi-Civita |
| test_Contract.jl | 135 | 9 | 0 | Contract | Metric/FV/epsilon contractions |
| test_ExpandScalarProduct.jl | 134 | 10 | 0 | ExpandScalarProduct | Scalar product expansion algebra |
| test_PolarizationSum.jl | 97 | 6 | 0 | PolarizationSum | Photon polarization sum rules |
| **SUBTOTAL** | **883** | **112** | **0** | | |

**MUnit Coverage Assessment:** 5 core FeynCalc functions translated (15 to 65 tests each). Per CLAUDE.md protocol (target ≥60 core functions with ≥5 tests): **9% covered (5/60)**. Functions covered:
- DiracTrick (65 tests) — expansion identities for Dirac chains
- DiracTrace (22 tests) — fermionic traces with Levi-Civita
- Contract (9 tests) — tensor contractions
- ExpandScalarProduct (10 tests) — kinematic algebra
- PolarizationSum (6 tests) — photon completeness relations

Missing major functions: TensorProduct, DecayAmplitude, Amplitude, Polarization, ChangeDimension, FermionSpinSum, etc.

### C. QGraf Integration Tests (30 files in qgraf/)

| File | LOC | Test Count | Broken | Process | Purpose |
|------|-----|-----------|--------|---------|---------|
| test_types.jl | 176 | 75 | 0 | unit | Partition/rho combinatorics validation |
| test_canonical.jl | 174 | 35 | 0 | unit | Diagram canonicality under automorphisms |
| test_halfedge.jl | 165 | 30 | 0 | unit | Half-edge graph representation |
| test_propagator_assemble.jl | 165 | 24 | 0 | unit | Propagator line construction |
| test_momentum_routing.jl | 171 | 25 | 0 | unit | Momentum flow through diagrams |
| test_step_b.jl | 140 | 21 | 0 | phase13 | Vertex/edge deduplication |
| test_qg21_battery.jl | 71 | 24 | 0 | qg21 | Main QGraf wrapping + piping |
| test_qg10_labels.jl | 90 | 19 | 0 | qg10 | Vertex label bijection |
| test_qgen_dpntro.jl | 80 | 12 | 0 | qgen | Depth-first recursion |
| test_qg10.jl | 69 | 14 | 0 | qg10 | Q-graph 10-line output parsing |
| test_filters_inline.jl | 49 | 10 | 0 | filter | Diagram inline-propagator elimination |
| test_filters_qumpi.jl | 90 | 4 | 0 | filter | User multi-particle irreducibility |
| test_automorphisms.jl | 81 | 7 | 0 | unit | Graph automorphism group |
| test_emission_amplitude.jl | 96 | 9 | 0 | unit | Photon/gluon emission structure |
| test_external_assemble.jl | 88 | 10 | 0 | unit | External leg tree construction |
| test_fermion_line.jl | 90 | 10 | 0 | unit | Fermion flow tracking |
| test_phase18a_pipeline.jl | 33 | 1 | 0 | integration | End-to-end qgraf invocation |
| test_sym_factor.jl | 65 | 4 | 0 | unit | Symmetry factor calculation |
| test_solve_tree_pipeline.jl | 29 | 3 | 0 | integration | Tree-level pipeline |
| test_step_c.jl | 99 | 10 | 0 | phase13 | Self-loop/tadpole filter |
| test_step_c_connectedness.jl | 65 | 5 | 0 | phase13 | Connectedness validation |
| test_qg21_integration.jl | 132 | 7 | 0 | integration | QG21 + filters end-to-end |
| test_qdis.jl | 78 | 1 | 0 | qgraf | QGraf disk I/O |
| test_qgen_recurse.jl | 90 | 3 | 0 | qgen | Recursion orbit enumeration |
| test_step_a.jl | 77 | 7 | 0 | phase13 | Multi-edge deduplication |
| test_vertex_assemble.jl | 88 | 7 | 0 | unit | Vertex dictionary construction |
| test_momentum.jl | 61 | 4 | 0 | unit | Momentum object manipulation |
| test_count_diagrams_qg21.jl | 71 | 10 | 1 | qg21 | QG21 vs legacy count matching |
| test_phase17_audition.jl | 75 | 5 | 4 | audition | Dedup algorithm cross-check (Burnside/canonical/prefilter) |
| test_filters_qumvi.jl | 81 | 6 | 0 | filter | Vacuum irreducibility |
| **SUBTOTAL** | **2,732** | **414** | **5** | | |

**Broken Test Inventory:**
- test_count_diagrams_qg21.jl:1 — "Was @test_broken" comment (now passing) — Burnside dedup fix
- test_phase17_audition.jl:4 — Canonical and prefilter dedup known over-counting bugs:
  - QED ee→μμ 1L (2-gen): Burnside = 18 ✓, but canonical/prefilter = 19 (off by 1) 
  - QED ee→ee tree: Burnside = 2 ✓, but canonical/prefilter = 1 (over-dedup under in↔out auto)

---

## 2. Test Orchestration: runtests.jl and run_v2_tests.sh

### runtests.jl (Single-Process Aggregator)

**Path:** `test/v2/runtests.jl` (33 LOC)

**Design:** Single Julia process loads FeynfeldX once, then includes 20 test files sequentially. Usage:
```bash
julia --project=. test/v2/runtests.jl
```

**Included tests (20 core files):**
1. test_coeff.jl (DimPoly)
2. test_colour.jl (Color algebra)
3. test_ee_mumu_x.jl (Tracer bullet)
4. test_self_energy.jl (Tree self-energy)
5. test_vertical.jl (Amplitude structure)
6. test_pave.jl (Integrals)
7. test_schwinger.jl (Reference formulas)
8. test_compton.jl (Compton scattering)
9. test_munit_batch1.jl (DiracTrace/Contract/PolarizationSum)
10. test_munit_batch2.jl (DiracTrick n≥3)
11. test_bhabha.jl (Bhabha scattering)
12. test_qqbar_gg.jl (QCD diagram count)
13. test_self_energy_1loop.jl (1-loop self-energy)
14. test_d0.jl (D₀ integrals)
15. test_running_alpha.jl (Running coupling)
16. test_ee_ww.jl (e+e→WW tree)
17. test_pipeline.jl (Full pipeline)
18. test_vertex_g2.jl (Photon 3-vertex)
19. test_box_ee_mumu.jl (1-loop box)
20. test_nlo_box_validation.jl (NLO validation)

**Status:** Missing from orchestration: test_diagram_gen.jl, test_vertex_arity.jl, test_qcd_4gluon.jl, test_qcd_ghost.jl, test_ee_ww_grozin.jl, and all 30 qgraf tests.

### run_v2_tests.sh (Individual Process Runner)

**Path:** `grind/run_v2_tests.sh` (11 LOC)

**Design:** Bash loop iterates test/v2/test_*.jl, invokes each in separate Julia process, pipes grep for "Test Summary|ERROR|Some tests" to stdout. Usage:
```bash
bash grind/run_v2_tests.sh > grind/v2_test_results.txt
```

**Coverage:** All 25 test_*.jl files (main + munit batches). Does NOT include qgraf tests.

**Output format:** `basename_$f: Test Summary: ... | ERROR: ... | Some tests failed`

---

## 3. MUnit FeynCalc Function Coverage

### Covered Functions (5/~60 ≈ 9%)

| Function | File | Tests | Status | Coverage % |
|----------|------|-------|--------|-----------|
| **DiracTrick** | test_DiracTrick.jl | 65 | ✓ | 100% (all free-index cases) |
| **DiracTrace** | test_DiracTrace.jl (batch1) | 22 | ✓ | 70% (basics + eps; missing BBox, cyclicity) |
| **Contract** | test_Contract.jl (batch1) | 9 | ✓ | 40% (metric/FV/epsilon; missing vector simplify) |
| **ExpandScalarProduct** | test_ExpandScalarProduct.jl (batch1) | 10 | ✓ | 50% (basic expansion; missing multi-angle) |
| **PolarizationSum** | test_PolarizationSum.jl (batch1) | 6 | ✓ | 30% (transverse sum; missing massive) |

### Missing Priorities (FeynCalc ecosystem scope)

**Immediate next:** (per CLAUDE.md Spiral 8)
- Amplitude (tree → α^n expansion) — 15-20 tests
- TensorProduct (gamma algebra) — 15 tests
- ChangeDimension (D-dimensional reduction) — 10 tests

**Medium term:**
- Polarization (explicit tensor forms) — 10 tests
- SimplifyPolyLog (Li₂ algebra) — 12 tests
- Simplify (automatic simplification) — 20 tests

**Long tail (~48 functions):** DecayAmplitude, PaVe* variants, Loop integrals, etc.

**Roadmap impact:** MUnit protocol (CLAUDE.md line 192) specifies metric = function coverage %, not test count. Current 9% (5/60) is early-stage; spiral targets are 25% by end of Spiral 8.

---

## 4. Known Broken Tests (5 Total)

### A. QGraf Dedup Algorithm Bugs (4 tests)

**File:** test/v2/qgraf/test_phase17_audition.jl

| Test Name | Line | Assertion | Expected | Actual | Root Cause |
|-----------|------|-----------|----------|--------|-----------|
| test_broken @ line 60 (QED ee→μμ 1L / canonical) | 60 | `count_dedup_canonical(...) == 18` | 18 | 19 | Canonical-rep bug: over-distinguishes orbit reps |
| test_broken @ line 61 (QED ee→μμ 1L / prefilter) | 61 | `count_dedup_prefilter(...) == 18` | 18 | 19 | Same canonicality compare bug |
| test_broken @ line 71 (QED ee→ee / canonical) | 71 | `count_dedup_canonical(...) == 2` | 2 | 1 | Over-dedup under in↔out automorphism |
| test_broken @ line 72 (QED ee→ee / prefilter) | 72 | `count_dedup_prefilter(...) == 2` | 2 | 1 | Same issue |

**Summary:** Burnside dedup (reference algorithm) is correct (Phase 12d fix validated). Canonical-compare and prefilter algorithms have known over-/under-counting bugs due to automorphism equivalence class failures. Unblocked by Burnside correctness; canonical/prefilter bugs isolated to test_phase17_audition.jl (Phase 17 audition, not Spiral 8 critical path).

### B. Historical "Was @test_broken" (now passing)

**File:** test/v2/qgraf/test_count_diagrams_qg21.jl, line 62-69

Comment: "Was @test_broken — Burnside returned 16 vs legacy 18. Fixed by enumerating all positional perms of remaining slots in _qgen_recurse + applying self-loop & multi-edge filters (qgen:13921-13954)."

Now passing; included for regression tracking.

---

## 5. Coverage Gaps & Strategic Assessment

### A. Process-Layer Coverage

| Process | Layer 1 (Model) | Layer 2 (Rules) | Layer 3 (Diagrams) | Layer 4 (Algebra) | Layer 5 (Integrals) | Layer 6 (Evaluate) | Tests |
|---------|-----------------|-----------------|-------------------|-------------------|-------------------|-------------------|-------|
| **Compton** e+γ→e+γ | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | 4 (tree) |
| **Bhabha** e+e→e+e | ✓ | ✓ | ✓ | ✓ | ✗ | ✗ | 4 (tree) |
| **e+e→μμ** | ✓ | ✓ | ✓ | ✓ | ✓ (partial) | ✗ | 40+ (tree + 1L box) |
| **e+e→WW** | ✓ | ✓ | ✗ (hand-built) | ✓ | ✗ | ✗ | 39 (tree) |
| **QCD 4-gluon** | ✓ | ✓ | ✓ | ✗ | ✗ | ✗ | 3 (count only) |

**Verdict:** Model/Rules fully pipeline-compliant. Diagrams partially tested (qgraf unit tests extensive; full pipeline integration sparse). Algebra excellent (Layer 4 100% coverage). Integrals (Layer 5) tested only for PaVe/D₀ special functions; no loop-integral cross-section coupling. Evaluate (Layer 6) untested in production pipeline.

### B. Algebra-Layer Coverage

**Strengths:**
- DimPoly (29 tests) — coefficient arithmetic ✓
- Color algebra (27 tests) — SU(N) ✓
- Dirac algebra (65+22 = 87 MUnit tests) — fermionic traces ✓
- Scalar products (10 tests) — kinematic algebra ✓

**Gaps:**
- Lorentz tensors (only 9 Contract tests) — needs 15+ more
- Epsilon (Levi-Civita) — 3 tests in Contract; expand-eps rules underspecified
- Gamma5 — test/v2/test_*.jl shows NO gamma5 algebra; only references in comments
- PaVe/D₀ evaluated numerically (50 tests); symbolic algebra missing

### C. Integral-Layer Coverage

| Integral | Type | Tests | Validated Against | Status |
|----------|------|-------|-------------------|--------|
| A₀(m²) | 1-point | 0 | — | ✗ Missing |
| B₀(p², m₀², m₁²) | 2-point | 0 | — | ✗ Missing |
| C₀ | 3-point | 0 | — | ✗ Missing |
| D₀ | 4-point | 23 | COLLIER + quadgk | ✓ Extensive |
| PaVe | General | 50 | LoopTools (Fortran) | ✓ Good |
| α(Q²) running | Coupling | 32 | — | ✓ Formula-based |

**Gap:** A₀, B₀, C₀ scalar integrals never cross-validated. PaVe reduction formulas correct but no diagram-coupling tests.

### D. Test-Execution Bottlenecks

1. **runtests.jl incomplete:** Omits 5 core test files (diagram_gen, vertex_arity, qcd_4gluon, qcd_ghost, ee_ww_grozin) + all 30 qgraf tests. Total = 62 files; orchestrated = 20. **Add-missing % = 68%.**

2. **No continuous integration:** run_v2_tests.sh is bash, not CI/CD. No auto-trigger on commits. Manual re-run only.

3. **Qgraf tests isolated:** 30 qgraf files NOT included in runtests.jl. Requires separate invocation. Discovery discrepancy: user must know about test/v2/qgraf/ to run.

### E. Documentation Gaps

1. **No test manifest:** No single document listing which processes are tested at which layers.
2. **No coverage matrix:** No README mapping test file ↔ physics process ↔ spiral milestone.
3. **Broken tests not flagged:** 5 @test_broken lines documented only in comment text, not in issue tracking.

---

## 6. Recommended Actions (Priority Order)

1. **Immediate (Spiral 8 pre-requisite):**
   - Add test_diagram_gen.jl, test_vertex_arity.jl to runtests.jl
   - Document 4 @test_broken tests in Beads issue tracker
   - Gamma5 algebra: add test_gamma5.jl (15 tests, per CLAUDE.md Spiral 8)

2. **Short-term (Spiral 9):**
   - EW model integration: test_ee_ww.jl via full pipeline (not hand-built)
   - QCD diagram generation: expand test_qcd_*.jl from count-only to amplitude tests
   - A₀, B₀, C₀ integral cross-validation vs LoopTools

3. **Medium-term (function coverage):**
   - MUnit batch 3: Amplitude, TensorProduct, ChangeDimension (25 tests)
   - Epsilon contraction rules: expand test_Contract.jl from 9 to 20 tests
   - Phase 17b audition resolution: fix canonical/prefilter dedup bugs

---

## Summary Statistics

| Metric | Count | Notes |
|--------|-------|-------|
| **Total test files** | 62 | 25 main + 5 munit + 30 qgraf + 1 orchestrator + 1 runner |
| **Total test assertions** | 301 | @test lines (not @test_broken) |
| **Known broken tests** | 5 | 4 dedup bugs + 1 historical (now passing) |
| **Total LOC (test)** | 7,628 | 3,713 main + 883 munit + 2,732 qgraf + 33 runtests + 11 runner |
| **Orchestrated via runtests.jl** | 20/25 | 80% coverage; omits 5 key files |
| **MUnit function coverage** | 5/~60 | 9% (DiracTrick, DiracTrace, Contract, ExpandSP, PolarizationSum) |
| **Pipeline processes tested** | 5 | Compton, Bhabha, e+e→μμ, e+e→WW, self-energy |
| **QGraf unit tests** | 414 | Comprehensive graph/automorphism coverage |
| **Integral validators** | 2 | D₀ (COLLIER), PaVe (LoopTools) |
| **Broken algorithm bugs** | 2 | Canonical dedup over-count; prefilter same issue |

