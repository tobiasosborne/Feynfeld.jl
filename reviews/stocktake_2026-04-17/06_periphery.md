# Periphery Stocktake — 2026-04-17

Covers everything outside the six stocktake "main acts" (`01_algebra`, `02_model_rules_diagrams`, `03_integrals_evaluate`, `04_qgraf_port`, `05_tests`): the frozen v1 codebase, top-level scripts, `grind/` diagnostic infrastructure, the `reviews/` scratch directory, the `refs/` reference corpus, and the top-level documentation files.

**Key takeaway:** v1 (`src/algebra`, `src/integrals`, 27 source files, 3,005 LOC; 20 test files, 2,344 LOC) is entirely superseded by v2. Every v1 file has a v2 counterpart except a handful of dead module scaffolds (`src/model/`, `src/rules/`, `src/diagrams/`, `src/evaluate/` and the matching `test/` subdirectories are all empty). `JULIA_PATTERNS.md` duplicates CLAUDE.md §6 verbatim (cosmetic formatting diffs only). `reviews/` is gitignored and holds 6 architectural review reports + 14 research/ids notes from the 2026-03-29 six-agent audit — useful as historical context but safe to archive.

---

## 1. v1 Source Inventory (`src/Feynfeld.jl` + `src/algebra/` + `src/integrals/`)

The v1 module entry `src/Feynfeld.jl` (75 LOC) includes all 23 `algebra/*.jl` and all 4 `integrals/*.jl` files; the Model/Rules/Diagrams/Evaluate includes are commented out.

| v1 file | LOC | v2 counterpart | Notes |
|---|---:|---|---|
| `src/Feynfeld.jl` | 75 | `src/v2/FeynfeldX.jl` (157) | v1 module entry; v2 adds Model/Rules/Diagrams/loop/phi3/QCD/EW includes. |
| `src/algebra/alg_expr.jl` | 166 | `src/v2/expr.jl` (150) | AlgTerm/AlgSum tree; v2 renamed. |
| `src/algebra/alg_ops.jl` | 194 | merged into v2 `contract.jl`, `expand_sp.jl`, `dirac_trace.jl` | v2 dispatches on AlgSum directly in each op file. |
| `src/algebra/colour_simplify.jl` | 116 | `src/v2/colour_simplify.jl` (183) | SU(N) simplifier; v2 grew for Spiral 6. |
| `src/algebra/colour_trace.jl` | 84 | `src/v2/colour_trace.jl` (84) | 1:1 port. |
| `src/algebra/colour_types.jl` | 148 | `src/v2/colour_types.jl` (122) | v2 trimmed. |
| `src/algebra/contract.jl` | 156 | `src/v2/contract.jl` (212) | v2 dispatch-native. |
| `src/algebra/dimensions.jl` | 84 | merged into v2 `coeff.jl` (DimPoly) | v2 replaced ad-hoc dim tracking with DimPoly. |
| `src/algebra/dirac_chain.jl` | 100 | v2 has no chain type — DOT integrated into `dirac_expr.jl` (107) | v2 dropped the wrapper. |
| `src/algebra/dirac_equation.jl` | 85 | *not in v2* (see ids_capabilities.txt: feynfeld-a0j) | v1 pattern kept as reference; v2 reimplementation pending. |
| `src/algebra/dirac_order.jl` | 65 | *not in v2* (see ids_capabilities.txt) | reference-only. |
| `src/algebra/dirac_scheme.jl` | 46 | *not in v2* (gamma5 scheme) | Spiral 8 scope. |
| `src/algebra/dirac_simplify.jl` | 44 | *not in v2* (orchestrator missing; feynfeld-isx) | reference-only. |
| `src/algebra/dirac_trace.jl` | 136 | `src/v2/dirac_trace.jl` (158) | v2 no-gamma5; else feature-parity. |
| `src/algebra/dirac_trick.jl` | 232 | `src/v2/dirac_trick.jl` (186) | v2 trimmed. |
| `src/algebra/dirac_types.jl` | 155 | `src/v2/dirac.jl` (120) | v2 parametric types. |
| `src/algebra/eps.jl` | 150 | `src/v2/eps_contract.jl` (86) + `src/v2/eps_evaluate.jl` (65) | v2 split contract/eval. |
| `src/algebra/expand_sp.jl` | 62 | `src/v2/expand_sp.jl` (95) | feature parity. |
| `src/algebra/fermion_spin_sum.jl` | 103 | `src/v2/spin_sum.jl` (171) + `src/v2/polarization_sum.jl` (57) | v2 split fermion/boson. |
| `src/algebra/minkowski.jl` | 47 | *not in v2* — TensorGR bridge removed | v2 no longer imports TensorGR. |
| `src/algebra/momentum.jl` | 137 | merged into `src/v2/types.jl` (101) | Momentum/MomentumSum in v2 types. |
| `src/algebra/pair.jl` | 137 | `src/v2/pair.jl` (75) | v2 leaner. |
| `src/algebra/sp_context.jl` | 54 | `src/v2/sp_context.jl` (71) + `src/v2/sp_lookup.jl` (32) | v2 split context/lookup. |
| `src/algebra/types.jl` | 82 | `src/v2/types.jl` (101) | merged with momentum.jl. |
| `src/integrals/feynamp_denominator.jl` | 122 | *not in v2* | v1 FAD not ported (v2 uses explicit propagators). |
| `src/integrals/pave.jl` | 86 | `src/v2/pave.jl` (86) | 1:1 port. |
| `src/integrals/pave_reduce.jl` | 56 | `src/v2/tid.jl` (175) | v2 renamed to TID (tensor-integral decomposition). |
| `src/integrals/tdec.jl` | 83 | `src/v2/tid.jl` (175) | merged into TID. |
| `src/model/`, `src/rules/`, `src/diagrams/`, `src/evaluate/` | 0 | v2 has `model.jl`, `rules.jl`, `diagrams.jl`, `diagram_gen.jl`, `cross_section.jl` etc. | **confirmed empty** — dead module scaffolds. |

**v1 total:** 27 files, **3,005 LOC**. No active includes of Model/Rules/Diagrams/Evaluate — those were never populated in v1.

---

## 2. v1 Test Inventory (`test/algebra/`, `test/integrals/`, `test/runtests.jl`)

The v1 entry `test/runtests.jl` (69 LOC) drives 18 algebra tests, 2 integrals tests, and 1 top-level `test_ee_mumu.jl`.

| v1 test file | LOC | v2 counterpart | Notes |
|---|---:|---|---|
| `test/runtests.jl` | 69 | `test/v2/runtests.jl` (33) | v2 entry orders 20 files. |
| `test/test_ee_mumu.jl` | 150 | `test/v2/test_ee_mumu_x.jl` (141) + `test/v2/test_pipeline.jl` (372) | v2 splits unit + pipeline. |
| `test/algebra/test_alg_expr.jl` | 233 | merged into `test/v2/test_coeff.jl` (64), `test_colour.jl` (144) etc. | v2 distributes. |
| `test/algebra/test_colour.jl` | 108 | `test/v2/test_colour.jl` (144) | expanded. |
| `test/algebra/test_contract.jl` | 160 | `test/v2/munit/test_Contract.jl` (MUnit-format) | v2 uses FeynCalc MUnit naming. |
| `test/algebra/test_dimensions.jl` | 52 | merged into `test/v2/test_coeff.jl` | v2 tests DimPoly. |
| `test/algebra/test_dirac_chain.jl` | 82 | *not in v2* — chain type dropped | obsolete. |
| `test/algebra/test_dirac_phase1b.jl` | 137 | covered by `test/v2/munit/test_DiracTrick.jl` + `test_DiracTrace.jl` | obsolete name. |
| `test/algebra/test_dirac_trace.jl` | 68 | `test/v2/munit/test_DiracTrace.jl` | MUnit port. |
| `test/algebra/test_dirac_trick.jl` | 98 | `test/v2/munit/test_DiracTrick.jl` | MUnit port. |
| `test/algebra/test_eps.jl` | 92 | covered by `test/v2/munit/test_Contract.jl` + research notes | partially ported. |
| `test/algebra/test_expand_sp.jl` | 148 | `test/v2/munit/test_ExpandScalarProduct.jl` | MUnit port. |
| `test/algebra/test_fermion_spin_sum.jl` | 120 | `test/v2/munit/test_PolarizationSum.jl` | MUnit port + expansion. |
| `test/algebra/test_minkowski.jl` | 50 | *not in v2* (TensorGR bridge removed) | obsolete. |
| `test/algebra/test_momentum.jl` | 35 | covered by `test/v2/qgraf/test_momentum.jl` (107) + `test_momentum_routing.jl` | v2 expanded. |
| `test/algebra/test_munit_dirac.jl` | 129 | `test/v2/test_munit_batch1.jl` (183) + `test/v2/test_munit_batch2.jl` (127) | v2 batch rename. |
| `test/algebra/test_munit_lorentz.jl` | 159 | split between `test/v2/test_munit_batch*.jl` and `munit/test_*.jl` | v2 reorganised. |
| `test/algebra/test_pair.jl` | 139 | covered incidentally in v2 pipeline tests | no dedicated v2 pair test. |
| `test/algebra/test_sp_context.jl` | 41 | covered incidentally | no dedicated v2 test. |
| `test/algebra/test_types.jl` | 123 | covered incidentally | no dedicated v2 test. |
| `test/integrals/test_pave.jl` | 103 | `test/v2/test_pave.jl` (184) | v2 expanded. |
| `test/integrals/test_tdec.jl` | 48 | `test/v2/test_pave.jl` (184) + `test/v2/test_d0.jl` (107) | merged. |
| `test/references/`, `test/rules/`, `test/model/`, `test/diagrams/`, `test/evaluate/` | 0 | v2 equivalents live in `test/v2/test_*.jl` and `test/v2/qgraf/` | **confirmed empty**. |

**v1 test total:** 21 files, **2,344 LOC**. All obsolete or already ported to v2.

---

## 3. Scripts Inventory (`scripts/`)

| Script | LOC (approx) | Purpose |
|---|---|---|
| `scripts/audition_compare.jl` | ~80 | Compares the three dedup strategies (Burnside / canonical / prefilter) across a battery of QED/φ³ tree+1L cases. Used during Phase 17 audition. |
| `scripts/debug_qed_1l.jl` | ~40 | Ad-hoc diagnostic: runs qg21 → qg10 → qgen stack for QED ee→μμ 1L and prints per-topology assignment counts. |
| `scripts/qgraf_golden_master_report.jl` | ~110 | Runs every supported golden-master case through legacy `count_diagrams`, prints PASS/FAIL/SKIP. Honours `QGRAF_MAX_SECONDS`. Beads `feynfeld-ney`. |
| `scripts/qgraf_golden_master_qg21_report.jl` | ~110 | Same pattern for Strategy C (`count_diagrams_qg21`). Sibling diagnostic. |

All four shell out to `src/v2/FeynfeldX.jl`. None touches v1.

---

## 4. `grind/` — qgraf diagnostic infrastructure

Per `grind/README.md`: setup used in HANDOFF Session 22 to find the qgen flavor-loop under-count bug. Workflow: copy qgraf's `qgraf-4.0.6.f08` into `grind/qgraf-grind.f08`, add `WRITE(0, ...)` instrumentation at qgen entry/match/emit points, build with `gfortran -O0 -g`, mirror in `dump_julia_emissions.jl`, diff the traces.

**Committed (our work, per `.gitignore` allow-list):**

| File | Purpose |
|---|---|
| `README.md` | workflow docs |
| `compare_topos.jl` | topology diff helper |
| `ctrl.dat`, `ctrl_phi3_2L.dat` | qgraf control files for ee→μμ 1L and φ³ 2L |
| `dump_julia_emissions.jl`, `dump_julia_phi3_2L.jl` | Julia diagnostic dumpers |
| `parse_grind_trace.jl` | parse qgraf WRITE(0,...) output |
| `inspect_dpntro.gdb`, `inspect_qgen.gdb` | gdb scripts |
| `run_v2_tests.sh` | sequential v2 test runner |

**Gitignored (regenerable):**

- `qgraf-grind`, `qgraf-grind.f08`, `fmodules/`, `models/`, `styles/`, `out_tmp/`, `qgraf-grind` binary — all derived from licensed qgraf source.
- `grind_trace.txt`, `grind_stdout.txt`, `julia_trace.txt`, `julia_stdout.txt`, `julia_phi3_2L.txt`, `grind_phi3_2L.txt`, `qgen_trace.txt`, `v2_test_results.txt` — diagnostic logs.

All diagnostics reproducible from the committed scripts + qgraf source.

---

## 5. `output/` — COLLIER evaluator output

- `ErrOut.cll`, `ErrOut.coli`, `InfOut.cll` — runtime logs from COLLIER's Fortran backend. Gitignored (see `.gitignore` line 55: `output/`). Regenerable by re-running D₀ evaluations. **Safe to delete.**

---

## 6. Reviews Inventory (`reviews/`, gitignored)

### Six architectural reviews (six-agent audit, 2026-03-29)

| File | One-line verdict |
|---|---|
| `01_architecture_review.md` | Core algebra (Layer 4) excellent; Layers 1-3 are scaffolding (bifurcation problem); Layer 5 has process-specific standalone recipes that bypass the pipeline. |
| `02_julia_idiomaticity_review.md` | "Very good foundation" — 6 CRITICAL + 9 HIGH items, mostly latent type instabilities (MomentumSum, gamma_pair, pair(), spin_sum Tuple{Any,...}, QEDModel.params Dict{Symbol,Any}, global mutable `_COLOUR_DUMMY_COUNTER`). |
| `03_vision_plan_review.md` | PRD vision realistic for QED/QCD/EW tree+1L scope; aspirational for SUSY/QG; sharpen MVP-for-external-users bar. |
| `04_reference_comparison.md` | ~41 of ~14,948 FeynCalc MUnit tests translated (0.27%); 5 dangerous gaps: gamma5 trace, Eps contraction, DiracEquation, D0 evaluation, diagram generation. |
| `05_test_quality_review.md` | 301 tests / 3,400 LOC = 0.6 test-to-source ratio (low); 3 known bugs with no regression tests (colour_trace i², B0 imaginary, vertex iε); error paths untested. |
| `06_ground_truth_audit.md` | CLAUDE.md Rule 2 audit: 5 CRITICAL uncited formulas, 1 verbatim equation mismatch, 8 HIGH vague citations, 12 MEDIUM P&S-only. All cited reference files exist on disk. |

### `ids_*.txt` (beads-issue seed lists, 5 files, ~45 lines total)

Seeds for beads issue creation, distilled from the reviews. Each line: `<bead-id>\t<SEVERITY>: <title>`.

- `ids_bugs.txt` (6) — bugs: DiracExpr.+ unbounded growth, global colour counter, B0 iε, etc.
- `ids_capabilities.txt` (13) — missing features: gamma5, DiracEquation, DiracSimplify, etc.
- `ids_citations.txt` (10) — citation fixes: B1 wrong Eq. number, missing A₀/B₀ closed-form citations, etc.
- `ids_tests.txt` (9) — regression tests needed + strengthenings.
- `ids_types.txt` (7) — type instabilities: MomentumSum, QEDModel.params, Tuple{Any,...}, etc.

### `research_*.txt` (9 files, ~3,000 lines total)

Deep-research notes produced by read-only research agents, feeding Spiral 8 and Spiral 9:

- `research_diagram_plan.txt` (400) — minimal 2→2 diagram generation plan (3 channels × N fields).
- `research_eps_contraction.txt` (646) — Levi-Civita contraction algorithm notes.
- `research_eps_feyncalc.txt` (309) — FeynCalc's Eps implementation analysis.
- `research_eps_subst.txt` (19) — confirms `substitute_index` does not handle Eps.
- `research_gamma5_codebase.txt` (227) — gamma5 uses in current codebase.
- `research_gamma5_formulas.txt` (546) — gamma5 trace formulas (BMHV, NDR, Larin, West).
- `research_gamma5_tests.txt` (232) — FeynCalc DiracTrace gamma5 MUnit cases.
- `research_julia_diagram_patterns.txt` (457) — Julia-native diagram-gen idioms.
- `research_momentum_sum.txt` (166) — MomentumSum caller analysis (feeds feynfeld-60n).

### `stocktake_2026-04-17/` (this stocktake)

Already contains `01_algebra.md`, `02_model_rules_diagrams.md`, `03_integrals_evaluate.md`, `04_qgraf_port.md`, `05_tests.md`. This file is `06_periphery.md`.

---

## 7. Reference Corpus (`refs/`, gitignored except README + qgraf scripts)

| Subdir | Size | Purpose |
|---|---:|---|
| `refs/COLLIER/` | 14M | One-loop scalar/tensor integral library (Fortran); v2 calls it via `d0_collier.jl`. |
| `refs/FeynArts/` | 16M | Mathematica diagram-generation reference; consulted for topology enumeration patterns. |
| `refs/FeynCalc/` | 48M | **Primary porting oracle.** 15,222 MUnit tests in `Tests/`. Ground truth for Dirac/Lorentz/colour algebra. |
| `refs/FeynRules/` | 3.8M | Lagrangian-to-Feynman-rules reference (Mathematica). |
| `refs/FormCalc/` | 88M | One-loop amplitude generator (Mathematica + FORM); consulted for PV reduction patterns. |
| `refs/LoopTools/` | 2.5M | One-loop numerics (Fortran); consulted for D₀ implementation. |
| `refs/qgraf/` | 19M | **qgraf 4.0.6** (Nogueira) diagram generator. Cleanroom `ALGORITHM.md` + 104 golden-master fixtures in `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/` (committed). Source `.f08` is gitignored (licence). |

All `refs/*` except `README.md`, `qgraf/README.md`, `qgraf/ALGORITHM.md`, qgraf golden-master scripts+fixtures, and `papers/README.md` + `papers/fetch_nogueira.mjs` are gitignored. Re-acquisition instructions in `refs/README.md`.

---

## 8. `refs/papers/` Inventory (all gitignored except README + fetch script)

| File | Reference | Citation role in source |
|---|---|---|
| `Peskin_Schroeder_1995.djvu` | Peskin & Schroeder, "An Introduction to QFT", Westview 1995 | Default textbook for Dirac/Lorentz/QED derivations. |
| `Nogueira1993_JCompPhys105_279.pdf` | Nogueira, J. Comput. Phys. 105 (1993) 279 | qgraf algorithm paper — cited in qgraf_port/ cleanroom. |
| `Denner1993_FortschrPhys41.pdf` | Denner, Fortschr. Phys. 41 (1993) 307 | One-loop techniques / PV reduction — cited in `src/v2/pave_eval.jl`, `tid.jl`. |
| `tHooftVeltman1979_NuclPhysB153.pdf` | 't Hooft & Veltman, Nucl. Phys. B153 (1979) 365 | Scalar integrals — cited in `src/v2/pave.jl`, `b0_eval.jl`. |
| `PassarinoVeltman1979_NuclPhysB160.pdf` | Passarino & Veltman, Nucl. Phys. B160 (1979) 151 | PV decomposition — cited in `src/v2/tid.jl`, `pave.jl`. |
| `MertigBohmDenner1991_FeynCalc_CPC64.pdf` | Mertig, Böhm, Denner, Comput. Phys. Commun. 64 (1991) 345 | FeynCalc origin paper. |
| `Shtabovenko2016_FeynCalc9_1601.01167.pdf` | Shtabovenko et al., arXiv:1601.01167 | FeynCalc 9 conventions. |
| `Shtabovenko2020_FeynCalc93_2001.04407.pdf` | Shtabovenko et al., arXiv:2001.04407 | FeynCalc 9.3 conventions. |
| `Shtabovenko2024_FeynCalc10_2312.14089.pdf` | Shtabovenko et al., arXiv:2312.14089 | **FeynCalc 10** — primary conventions reference. |
| `vanOldenborgh1990_ZPhysC46.pdf` | van Oldenborgh, Z. Phys. C46 (1990) 425 | LoopTools origin — cited in `src/v2/d0_collier.jl`. |
| `PDG2024_list_electron.pdf`, `_muon`, `_tau`, `_z_boson` | Particle Data Group 2024 particle listings | Mass/width values. |
| `PDG2024_rev_phys_constants.pdf`, `_quark_masses`, `_standard_model` | PDG 2024 review chapters | Fundamental constants + SM parameters — cited in `src/v2/ew_parameters.jl`. |
| `PDG2024_sum_leptons.pdf`, `_quarks` | PDG 2024 summary tables | Lepton/quark mass summaries. |
| `fetch_nogueira.mjs` | Our Playwright fetch script | Committed. Re-downloads Nogueira (1993) via TIB VPN + authenticated browser session. |
| `README.md` | Committed | Documents fetch pattern. |

Acquisition instructions in `refs/papers/README.md`. Nogueira script is the only automated fetch; others are manual TIB-VPN downloads.

---

## 9. Top-level Docs Inventory

| File | LOC | Purpose |
|---|---:|---|
| `CLAUDE.md` | ~200 | **Authoritative project instructions.** Architecture, Tobias's 12 rules, MUnit protocol, Julia idiom cheatsheet (§6). Auto-loaded into every agent session. |
| `Feynfeld_PRD.md` | (not measured here) | **Vision + spiral plan.** Six-layer architecture; replaces Mathematica ecosystem; §6 also contains Julia cheatsheet (three-way replica). |
| `JULIA_PATTERNS.md` | 88 | **Duplicates CLAUDE.md §6** (the DO / DO NOT cheatsheet). Confirmed: `diff` shows only cosmetic whitespace differences (`struct X; field::T; end` collapsed vs. multiline block form). File header explicitly states "This is replicated in `Feynfeld_PRD.md` §6 and `CLAUDE.md`. All three copies must stay in sync." |
| `SPIRAL_9_PLAN.md` | 195 | **Spiral 9 design doc (2026-03-29).** Goal: every existing process flows through the full six-layer pipeline, no more hand-built amplitudes. Key insight: tree-level 2→2 has exactly 3 channel topologies (s, t, u) — diagram generation is filtering, not enumeration. Scope: ~200 LOC across 3 files (`channels.jl`, filtering hook, audit). **NOT** a FeynArts port. |
| `AGENTS.md` | 150 | **Codex-agent-facing quick reference.** Tells `codex`-family agents to use `bd` (beads) for issue tracking and warns about interactive shell aliases (`cp -f`, `mv -f`, `rm -f`). Includes BEADS integration block (v:1, profile:full). Complements CLAUDE.md (which is Claude-Code-facing). |
| `HANDOFF.md` | (not measured) | Rolling session log; Session 26–27 summarise current Phase 18 state. |
| `Project.toml` / `Manifest.toml` | — | Julia package metadata. |
| `LICENSE` | — | Project licence. |

---

## 10. Recommendations

### Safe to delete now

1. **`src/algebra/*.jl` (23 files, 2,930 LOC).** All 23 files have v2 counterparts (see §1 table). Blocker: `src/Feynfeld.jl` module still `include`s them — delete that file too.
2. **`src/integrals/*.jl` (4 files, 347 LOC).** Same story — all ported to v2 (`pave.jl`, `tid.jl`).
3. **`src/Feynfeld.jl` (75 LOC).** v1 module entry. Replace with a thin re-export of `FeynfeldX` if backward compat is wanted, else delete.
4. **Empty scaffolds `src/model/`, `src/rules/`, `src/diagrams/`, `src/evaluate/`.** Zero files. Delete directories.
5. **`test/runtests.jl` + `test/test_ee_mumu.jl` + `test/algebra/` (18 files) + `test/integrals/` (2 files).** 2,344 LOC of v1 tests, all with v2 counterparts (see §2). Re-point `Pkg.test()` at `test/v2/runtests.jl`.
6. **Empty test scaffolds `test/model/`, `test/rules/`, `test/diagrams/`, `test/evaluate/`, `test/references/`.** Zero files.
7. **`JULIA_PATTERNS.md`.** Duplicates CLAUDE.md §6 verbatim (modulo whitespace). Either delete and rely on CLAUDE.md, or keep as standalone agent-onboarding doc but mark it as the canonical copy and remove from CLAUDE.md / PRD (violates DRY).
8. **`output/*.cll`, `*.coli`.** Gitignored runtime logs; already excluded from git. Can be wiped any time.
9. **`grind/*_trace.txt`, `grind/*_stdout.txt`, `grind/v2_test_results.txt`.** Gitignored diagnostic logs. Regenerable.

### Keep

1. **v2 code and tests** (`src/v2/*.jl` 54 files, `test/v2/*.jl` + `test/v2/munit/` + `test/v2/qgraf/`). Active.
2. **`scripts/*.jl` (4 files).** All target v2; all currently useful (golden-master audit, Phase 17 audition compare, 1L QED debug).
3. **`grind/` committed files** (README, dump scripts, .dat files, .gdb scripts, `run_v2_tests.sh`, `parse_grind_trace.jl`, `compare_topos.jl`). The infrastructure for future qgraf-vs-Julia bisections. Cost: ~10 files, mostly small.
4. **`reviews/NN_*.md` (6 files) + `ids_*.txt` (5) + `research_*.txt` (9).** Historical record of the 2026-03-29 six-agent audit. Tied to specific beads issues (`feynfeld-ccy`, `-qde`, `-1rb`, `-d2g`, etc.). Archive (e.g. move under `reviews/archive_2026-03-29/`) rather than delete — the citations and research notes are valuable context for fixing the bugs they flagged.
5. **`reviews/stocktake_2026-04-17/`.** This stocktake.
6. **`refs/` tree.** Reference corpus — all ground truth. Gitignored so cost is local-disk only.
7. **`refs/papers/`.** All cited; Rule 2 requires local copies. Keep.
8. **Top-level docs `CLAUDE.md`, `Feynfeld_PRD.md`, `HANDOFF.md`, `SPIRAL_9_PLAN.md`, `AGENTS.md`, `LICENSE`, `Project.toml`, `Manifest.toml`.** All active.

### Migration sequence (suggested)

```
1. Confirm test/v2/runtests.jl covers every useful v1 test (spot-check).
2. Move reviews/NN_*.md → reviews/archive_2026-03-29/.
3. git rm src/Feynfeld.jl src/algebra/ src/integrals/ src/model/ src/rules/ src/diagrams/ src/evaluate/.
4. git rm test/runtests.jl test/test_ee_mumu.jl test/algebra/ test/integrals/ test/model/ test/rules/ test/diagrams/ test/evaluate/ test/references/.
5. Point Project.toml [test] at test/v2/runtests.jl; promote src/v2/FeynfeldX.jl → src/Feynfeld.jl (or add re-export stub).
6. Collapse JULIA_PATTERNS.md duplication: keep it as canonical, shorten CLAUDE.md §6 and PRD §6 to a pointer.
7. Re-run full v2 suite to confirm no hidden v1 dependency.
```

Estimated deletion: **~3,000 LOC v1 source + ~2,344 LOC v1 tests + ~88 LOC JULIA_PATTERNS.md duplication = ~5,400 LOC removed**, zero capability lost. The repo goes from "v1 FROZEN / v2 ACTIVE" bifurcation to a single-tree codebase matching the CLAUDE.md vision ("v1 will be deleted. Do not extend.").
