# HANDOFF — 2026-03-29 (End of Session 9, Spiral 9 complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Read `Feynfeld_PRD.md` — vision, endgame interaction, spiral plan, hardness scale
3. Read `src/v2/DESIGN.md` — type system, anti-patterns, Session 8 review findings
4. Read `SPIRAL_9_PLAN.md` — the concrete plan for pipeline consolidation
5. Run `bd ready` to see available work (if beads errors: `bd init --force --prefix feynfeld && bd backup restore`)
6. Run `for f in test/v2/test_*.jl; do julia --project=. "$f"; done` to verify 307 tests pass
7. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## THE PIPELINE PRINCIPLE

**The pipeline IS the architecture. If a process bypasses it, the architecture is incomplete.**

Every process MUST flow through: Model → Rules → Diagrams → Algebra → Integrals → Evaluate.
No hand-built amplitudes in test files. No standalone recipes that skip the pipeline.
Standalone formulas (schwinger.jl, ew_cross_section.jl, etc.) are REFERENCE IMPLEMENTATIONS
for cross-validation — not substitutes for the pipeline.

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **CITE EVERYTHING (tiered).**
   - New formula: full triple (local path + eq number + verbatim).
   - Routine MUnit permutation: source file + test ID sufficient.
   ```julia
   # Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.21)
   # "A₀(m) = m²(Δ - ln(m²/μ²) + 1)"
   ```
3. **ACQUIRE PAPERS BEFORE CODE.** arXiv → TIB VPN → playwright-cli.
4. **JULIA IDIOMATIC.** Dispatch, not isa cascades. No Any. No global mutable state.
5. **WORKFLOW (TIERED):**
   - <5 LOC: direct fix, no subagents.
   - <20 LOC: 1 research + 1 review.
   - >20 LOC: 3 research + 1 review.
6. **REVIEW.** Rigorous reviewer after every core change.
7. **TESTING.** Known bugs MUST have `@test_broken`. Targeted or full suite in background.
8. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.
9. **LOC LIMIT ~200.** No source file exceeds ~200 lines.
10. **NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT OVERVIEW

### What is Feynfeld.jl?

Julia-native, agent-facing, full-stack physics computation suite. Lagrangian →
cross-section in one `using Feynfeld`. Replaces FeynRules + FeynArts + FeynCalc +
FormCalc + LoopTools. See PRD for the full vision.

### Branch and code location

- **Branch:** `experimental/rebuild-v2`
- **v2 source:** `src/v2/` (33 files, ~4,000 LOC)
- **v2 tests:** `test/v2/` (17 files, 324 tests)
- **v1:** `src/algebra/`, `src/integrals/` — FROZEN, will be deleted. Do NOT extend.

---

## WHAT EXISTS (v2, 307 tests)

### Six-layer pipeline

```
Layer 1: Model      → qed_model()                    → QEDModel
Layer 2: Rules      → feynman_rules(model)            → FeynmanRules (callable)
Layer 3: Diagrams   → tree_diagrams(model, ...)       → Vector{FeynmanDiagram}
Layer 4: Algebra    → trace → contract → expand → eval → AlgSum (scalar)
Layer 5: Integrals  → PaVe{N}, evaluate(::PaVe; mu2)  → ComplexF64
Layer 6: Evaluate   → solve_tree(prob) → σ             → Float64
```

**WARNING:** Only e+e-→μ+μ- uses the full pipeline. All other processes bypass
Layers 1-3 and construct amplitudes by hand in test files. This is the #1
architectural problem. See `SPIRAL_9_PLAN.md` for the fix.

### Layer 4: Algebra (the core — excellent, validated)

| File | What | Key types/functions |
|------|------|-------------------|
| `coeff.jl` | DimPoly coefficient algebra | `DimPoly`, `DIM`, `Coeff`, `evaluate_dim` |
| `types.jl` | Physics indices, momenta | `LorentzIndex`, `Momentum`, `MomentumSum`, `momentum_sum()` |
| `colour_types.jl` | SU(N) types | `AdjointIndex`, `FundIndex`, `SUNT`, `SUNF`, `SUND`, deltas |
| `pair.jl` | Parametric Pair{A,B} | `MetricTensor`, `FourVector`, `ScalarProduct`, `pair()` |
| `expr.jl` | Dict-based AlgSum | `AlgSum`, `AlgFactor` (6-type union), `FactorKey`, `alg()` |
| `sp_context.jl` | Scalar product context | `SPContext`, `ScopedValues`, `evaluate_sp` |
| `contract.jl` | Lorentz contraction | `contract()`, `substitute_index()`, Eps index handling |
| `eps_contract.jl` | Levi-Civita contraction | `eps_contract()`, ε·ε = -det[pair(aᵢ,bⱼ)] |
| `expand_sp.jl` | MomentumSum bilinear expansion | `expand_scalar_product()` |
| `dirac.jl` | Dirac gamma types | `DiracGamma{S}`, `Spinor{K}`, `GA/GAD/GS/GA5/GA6/GA7` |
| `dirac_trace.jl` | Trace → AlgSum | `dirac_trace()` — handles γ5, projectors, arbitrary chains |
| `dirac_expr.jl` | Matrix-valued expressions | `DiracExpr`, simplify-on-add |
| `dirac_trick.jl` | γ^μ...γ_μ contraction | `dirac_trick()` (n=0,1,2,3,4 + general n≥5) |
| `spin_sum.jl` | Fermion spin sums | `spin_sum_amplitude_squared()` |
| `colour_trace.jl` | SU(N) trace | `colour_trace()` (concrete N, recursive) |
| `colour_simplify.jl` | δ contraction, f·f, d·d | `contract_colour()` via `_self_contract` dispatch |
| `polarization_sum.jl` | Photon/gluon/massive pol sum | `polarization_sum()`, `polarization_sum_massive()` |

### Layers 1-3, 5-6

| File | What |
|------|------|
| `model.jl` | `QEDModel`, `Field{Species}`, `GaugeGroup`, traits |
| `rules.jl` | `FeynmanRules` callable, vertex/propagator dispatch on species |
| `diagrams.jl` | Hard-coded e+e-→μ+μ- s-channel topology (SCAFFOLDING) |
| `pave.jl` | `PaVe{N}` parametric type, named constructors A0/B0/B1/C0/C1/C2/D0 |
| `pave_eval.jl` | `evaluate(::PaVe{1,2,3})`, B0/C0 via QuadGK, C1/C2 via PV reduction |
| `vertex.jl` | REFERENCE: `vertex_f2_zero`, `vertex_f2` — QED g-2 (standalone) |
| `schwinger.jl` | REFERENCE: Schwinger correction formula + vacuum polarization |
| `running_alpha.jl` | REFERENCE: running α(q²), SM fermion table, delta_alpha |
| `ew_parameters.jl` | EW constants: M_W, M_Z, sin²θ_W, Z couplings (PDG 2024) |
| `ew_cross_section.jl` | REFERENCE: σ(e+e-→W+W-) Grozin formula (standalone) |
| `cross_section.jl` | `CrossSectionProblem`, `solve_tree`, Mandelstam |

### Test files (307 tests across 16 files)

| File | Tests | What |
|------|-------|------|
| `test_coeff.jl` | 29 | DimPoly arithmetic |
| `test_colour.jl` | 25 | SU(N) traces, δ contraction, f·f/d·d/f·d exact values |
| `test_ee_mumu_x.jl` | 14 | e+e-→μ+μ- algebra (P&S 5.10) |
| `test_self_energy.jl` | 25 | DiracExpr, DiracTrick n=0,1,2 |
| `test_vertical.jl` | 33 | Full pipeline: Model→Rules→Diagrams→Algebra→σ |
| `test_pave.jl` | 51 | PaVe types + A₀/B₀/B₁ numerical evaluation |
| `test_schwinger.jl` | 15 | Schwinger correction + vacuum polarization |
| `test_compton.jl` | 4 | Compton |M|² from pipeline vs P&S Eq. 5.87 |
| `test_munit_batch1.jl` | 23 | MUnit: DiracTrace, Contract, PolarizationSum |
| `test_munit_batch2.jl` | 18 | MUnit: DiracTrick n=3,4 |
| `test_bhabha.jl` | 4 | Bhabha |M̄|² at 2 kinematic points vs FeynCalc |
| `test_qqbar_gg.jl` | 2 | QCD qq̄→gg |M̄|² at 2 kinematic points |
| `test_self_energy_1loop.jl` | 13 | 1-loop Σ(p) via B₀, B₁, A₀ |
| `test_vertex_g2.jl` | 32 | C₀/C₁/C₂ evaluation, F₂(0)=α/(2π) |
| `test_running_alpha.jl` | 34 | Running α(q²), Δα, improved Born σ |
| `test_ee_ww.jl` | 36 | Tree-level e⁺e⁻→W⁺W⁻, massive pol sum exact values |

---

## WHAT WAS DONE IN SESSION 9

### Spiral 9 complete: Pipeline consolidation

**6 beads issues closed** (feynfeld-e8d, -czy, -2yw, -cdf, -07g + previous).

#### 1. Pipeline infrastructure (Phases A-B)

| File | LOC | What |
|------|-----|------|
| `channels.jl` (NEW) | 103 | `TreeChannel`, `tree_channels()` — enumerate s/t/u channels by vertex filtering |
| `amplitude.jl` (NEW) | 141 | `build_amplitude()` — boson exchange (2 chains) + fermion exchange (propagator decomp) |
| `diagrams.jl` (TRIMMED) | 22 | `ExternalLeg` only — old FeynmanDiagram/tree_diagrams removed |
| `cross_section.jl` (UPDATED) | — | `solve_tree` now uses tree_channels + build_amplitude |
| `spin_sum.jl` (EXTENDED) | +93 | `spin_sum_interference()` for reconnected traces, `_cross_line_trace()` |

#### 2. Model files (Phases C-D)

| File | LOC | What |
|------|-----|------|
| `qcd_model.jl` (NEW) | 60 | `QCDModel` with qqg + ggg vertices, `triple_gauge_vertex()` |
| `ew_model.jl` (NEW) | 49 | `EWModel` with eeγ/eeZ/eνW/WWγ/WWZ vertices |

#### 3. Pipeline tests

| File | Tests | What |
|------|-------|------|
| `test_pipeline.jl` (NEW) | 17 | Bhabha (3), Compton (4), qq→gg (2), ee→WW channels (8) |
| `test_vertical.jl` (UPDATED) | 34 | Uses tree_channels + build_amplitude (no more FeynmanDiagram) |

#### 4. Processes now through the pipeline

| Process | Channels | Exchange type | Status |
|---------|----------|---------------|--------|
| e+e-→μ+μ- | s(γ) | boson | Full pipeline (solve_tree) |
| Bhabha | s(γ)+t(γ) | boson+interference | Pipeline + manual compose |
| Compton | s(e)+u(e) | fermion | Pipeline + manual compose |
| qq̄→gg | t(q)+u(q)+s(g) | fermion+triple gauge | Pipeline for t/u, manual s |
| ee→WW | s(γ)+s(Z)+t(ν)+u(ν) | channel enumeration only | Model + channels only |

**Before Spiral 9**: Only e+e-→μ+μ- used the pipeline.
**After Spiral 9**: 4 of 5 processes have pipeline-generated amplitudes.

#### 5. Key design decisions

- `ExternalLeg.mass` field added (backward-compatible default 0//1)
- All amplitude indices use `DimD()` (D-dimensional, for consistent traces)
- `evaluate_m_squared` now calls `evaluate_dim` for DimPoly coefficients
- Channel-specific Lorentz indices (:mu_s, :mu_t, :mu_u) prevent interference collisions

---

## WHAT WAS DONE IN SESSION 8

### Spiral 8 complete: Architectural review + bug fixes + gamma5 + Eps

**24 issues closed** (79 → 103 of 145 total).

#### 1. Six-agent architectural review (reports in `reviews/`)

| Report | Focus | Key finding |
|--------|-------|-------------|
| `01_architecture_review.md` | Layer separation, type system | Core algebra excellent; Layers 1-3 are scaffolding |
| `02_julia_idiomaticity_review.md` | Antipatterns, code smells | 6 CRITICAL type instabilities found and fixed |
| `03_vision_plan_review.md` | Spiral methodology, risks | MUnit mop-up wrong move; diagram gen is #1 gap |
| `04_reference_comparison.md` | FeynCalc coverage map | 15% Dirac, 20% Lorentz, 25% colour coverage |
| `05_test_quality_review.md` | Coverage matrix, false positives | 7 MINIMAL-coverage files, several weak tests |
| `06_ground_truth_audit.md` | Citation verification | 6 CRITICAL uncited formulas, 1 wrong eq number |

#### 2. Type system fixes (all CRITICAL)

| Fix | File | What |
|-----|------|------|
| MomentumSum constructor | `types.jl` | Split into pure `new()` + `momentum_sum()` factory |
| gamma_pair return type | `dirac.jl` | Now returns `AlgSum` uniformly (eliminated 8 isa checks) |
| pair() factory | `pair.jl` | Split into two dispatch methods (LI×LI + general) |
| QEDModel.params | `model.jl` | Removed unused `Dict{Symbol, Any}` field |
| Tuple{Any,...} | `spin_sum.jl`, `expand_sp.jl` | Replaced with concrete types |
| _COLOUR_DUMMY_COUNTER | `colour_types.jl` | Replaced global Ref(0) with `gensym(:c)` |
| DiracExpr.+ | `dirac_expr.jl` | Now simplifies on add (prevents O(N²) blowup) |

#### 3. New capabilities

| Capability | File | LOC | What |
|-----------|------|-----|------|
| Eps contraction | `eps_contract.jl` (NEW) | 85 | ε·ε = -det[pair(aᵢ,bⱼ)] via Leibniz formula, 24 S₄ permutations |
| Eps index handling | `contract.jl` | +15 | `_indices(::Eps)`, `_subst_factor(::Eps)`, eps_contract pre-pass |
| gamma5 traces | `dirac_trace.jl` | +100 | Tr[γ5 γ^a γ^b γ^c γ^d] = -4ε^{abcd}, projector expansion, recursion |
| isa → dispatch | `colour_simplify.jl` | refactor | `_self_contract` dispatch methods, double-call fix |

#### 4. Citation audit (6 CRITICAL + 8 HIGH + 1 MEDIUM fixed)

All physics formulas now cite: local file path + equation number + verbatim.
Key fix: B₁ cited as "Denner Eq. B.7" was WRONG — correct is Eq. (B.9).

#### 5. Test strengthening (+6 tests)

- f·f = 3δ, d·d = (5/3)δ, f·d = 0 (exact value assertions)
- Massive pol sum: metric coeff = -1, k^μk^ν/M² coeff = 1/M² (exact)

#### 6. Architecture restructuring

- **Deferred 9 FeynArts/FeynRules-port subtasks** (described porting 3,500+ LOC Mathematica)
- **Created 6 new Julia-native architecture issues** with ~200 LOC total scope
- **Wrote `SPIRAL_9_PLAN.md`** — 3 files, 4 phases, concrete types
- **Updated CLAUDE.md** with pipeline principle, tiered workflow, tiered citations
- **Updated PRD** with endgame vision (pp→H+jet at NLO), hardness scale, revised spirals
- **Updated DESIGN.md** with Session 8 cockroaches and planned fixes

---

## WHAT TO DO NEXT

### Priority 1: Spiral 9 — Pipeline consolidation

The gateway is open. gamma5 traces and Eps contraction are implemented.
Read `SPIRAL_9_PLAN.md` for the full plan. Summary:

**Phase A: QED boson exchange** (~60 LOC)
- `feynfeld-e8d`: Channel type + `tree_channels()` — enumerate s/t/u channels
- `feynfeld-czy`: `build_amplitude()` — turn channels into DiracChains (blocked by e8d)
- File: `src/v2/channels.jl` (new)
- File: `src/v2/amplitude.jl` (new)
- Validation: ee→μμ through pipeline matches test_vertical.jl

**Phase B: QED fermion exchange** (~40 LOC)
- Extend `build_amplitude` for fermion propagator (p̸+m inside chain)
- Covers Compton + Bhabha t-channel
- `feynfeld-07g`: Pipeline test for ALL completed processes (blocked by czy)

**Phase C: QCD** (~40 LOC)
- `feynfeld-2yw`: QCDModel with qqg + ggg vertices
- Covers qq̄→gg through the pipeline

**Phase D: EW** (~60 LOC)
- `feynfeld-cdf`: EWModel with SM vertex dispatch (NOW UNBLOCKED)
- eeγ, eeZ (needs γ5 — done!), eνW, WWγ, WWZ vertices
- Covers ee→WW through the pipeline

### Priority 2: Remaining Spiral 8 bugs

| Issue | What | Priority |
|-------|------|----------|
| `feynfeld-83m` | colour_trace n≥4 missing i² factor | P1 (known bug, workaround exists) |
| `feynfeld-1rb` | B₀ imaginary part (iε prescription) | P1 (blocks timelike NLO) |
| `feynfeld-0h8` | @test_broken for colour_trace n≥4 | P1 (blocked by 83m) |

### Priority 3: Infrastructure

| Issue | What |
|-------|------|
| `feynfeld-akr` | CI/CD GitHub Actions |
| `feynfeld-qyu` | Rename FeynfeldX → Feynfeld |
| `feynfeld-zne` | Register package (blocked by akr + qyu) |

---

## KNOWN BUGS

- **`feynfeld-83m`** — colour_trace n≥4 missing i² factor. The recursive
  decomposition T^aT^b = δ^{ab}/(2N) + (1/2)(d+if)T^c omits the imaginary
  unit i on the f term. Workaround: use analytical colour factors (Casimirs).

- **B₀ imaginary part** (`feynfeld-1rb`) — `_B0_quadgk` returns approximate zero
  imaginary part for timelike momenta via `complex(f, -1e-30)` hack. Real part
  is correct to <1e-8. Needs proper analytic continuation.

- **vertex_f2 above threshold** — doesn't implement iε for q² > 4m². Use only
  for spacelike or below threshold.

---

## REFERENCE CODEBASES

All in `refs/` (gitignored):
- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`. Examples in `Examples/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — Loop integral numerics (Fortran source).
- `refs/papers/` — 16+ local paper copies. See `reviews/06_ground_truth_audit.md` for full list.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run all v2 tests (307 tests across 16 files)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Run specific test
julia --project=. test/v2/test_ee_ww.jl

# Fast smoke test (~2 min)
for f in test/v2/test_coeff.jl test/v2/test_colour.jl test/v2/test_vertical.jl \
         test/v2/test_schwinger.jl test/v2/test_running_alpha.jl test/v2/test_ee_ww.jl; do
    julia --project=. "$f"
done

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues
bd blocked            # dependency chains

# Session end protocol
git add <files> && git commit -m "..." && git push
bd backup save && bd dolt push
```

---

## FILE MAP

```
src/v2/
├── FeynfeldX.jl          # Module root, includes + exports
├── coeff.jl              # DimPoly coefficients
├── types.jl              # PhysicsIndex, Momentum, MomentumSum, momentum_sum()
├── colour_types.jl       # SU(N): AdjointIndex, FundIndex, SUNT, SUNF, SUND
├── pair.jl               # Parametric Pair{A,B}, dispatch-based pair() factory
├── expr.jl               # AlgSum (Dict), AlgFactor (6-type union), FactorKey
├── sp_context.jl         # SPContext + ScopedValues
├── contract.jl           # Lorentz contraction + Eps index handling + eps_contract pre-pass
├── eps_contract.jl       # NEW: ε·ε → -det[pair(aᵢ,bⱼ)] via Leibniz formula
├── expand_sp.jl          # Scalar product bilinear expansion
├── dirac.jl              # DiracGamma{S}, Spinor{K}, DiracChain, GA5/GA6/GA7
├── dirac_trace.jl        # Dirac trace → AlgSum (NOW WITH gamma5 + projectors)
├── spin_sum.jl           # Fermion spin sums
├── dirac_expr.jl         # DiracExpr: simplify-on-add (fixed Session 8)
├── dirac_trick.jl        # D-dim γ^μ...γ_μ for n=0..5+
├── colour_trace.jl       # SU(N) trace → AlgSum (BUG: i² missing for n≥4)
├── colour_simplify.jl    # Delta contraction via _self_contract dispatch
├── polarization_sum.jl   # Feynman + axial + massive gauge pol sums
├── model.jl              # AbstractModel, QEDModel, Field{Species}
├── rules.jl              # FeynmanRules callable
├── diagrams.jl           # FeynmanDiagram (SCAFFOLDING — to be replaced by channels.jl)
├── pave.jl               # PaVe{N} type, named constructors
├── pave_eval.jl          # evaluate(::PaVe) via QuadGK + PV reduction
├── schwinger.jl          # REFERENCE: Schwinger correction + vacuum pol
├── vertex.jl             # REFERENCE: QED vertex F₂(0)=α/(2π)
├── running_alpha.jl      # REFERENCE: SM fermions, running α(q²)
├── ew_parameters.jl      # EW SM: M_W, M_Z, sin²θ_W, Z couplings
├── ew_cross_section.jl   # REFERENCE: σ(e+e-→W+W-) Grozin formula
├── cross_section.jl      # Mandelstam, Problem/Solve, σ
├── DESIGN.md             # Design choices + Session 8 review findings
└── VERTICAL_PLAN.md      # Historical

test/v2/                  # 16 files, 307 tests
reviews/                  # 6 review reports + research files from Session 8
SPIRAL_9_PLAN.md          # Concrete plan: 3 files, ~200 LOC, 4 phases

Total: ~3,600 source LOC, ~2,100 test LOC, 307 tests, all files < 200 LOC.
```
