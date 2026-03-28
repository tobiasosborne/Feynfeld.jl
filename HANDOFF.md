# HANDOFF — 2026-03-28 (End of Session 5, Spiral 1 complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, anti-hallucination citation pattern, Julia idiom cheatsheet
2. Read `Feynfeld_PRD.md` — vision, spiral methodology, MUnit coverage target
3. Read `src/v2/DESIGN.md` — v2 type system, anti-patterns, cockroaches found
4. Read `JULIA_PATTERNS.md` — Julia idiom cheatsheet (same content as in CLAUDE.md §6)
5. Run `bd ready` to see available work (if beads errors, run `bd init --force --prefix feynfeld && bd backup restore`)
6. Run `for f in test/v2/test_*.jl; do julia --project=. "$f"; done` to verify 214 tests pass

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **ANTI-HALLUCINATION: CITE EVERYTHING.** Every formula must cite: local file
   path + equation number + verbatim copy. No exceptions. STOP if source missing.
   ```julia
   # Ref: refs/papers/Denner1993.pdf, Eq. (4.18)
   # "B₁(p², m₀², m₁²) = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)"
   ```
3. **JULIA IDIOMATIC ALL THE WAY.** Read the cheatsheet in CLAUDE.md FIRST.
   Use existing packages (QuadGK, PolyLog). No hand-rolled numerics. No isa cascades.
   Parametric types + dispatch. No OOP struct hierarchies unless dispatch genuinely needed.
4. **WORKFLOW.** 3 subagents before core code change: research source + 2 solutions.
5. **REVIEW.** Rigorous reviewer agent after every core change.
6. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.

**NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT OVERVIEW

### What is Feynfeld.jl?

Julia-native, agent-facing, full-stack physics computation suite. Lagrangian →
cross-section in one `using Feynfeld`. Replaces FeynRules + FeynArts + FeynCalc +
FormCalc + LoopTools. See PRD for the full vision (eventually all of physics:
SUSY, quantum gravity, etc.).

### Development methodology: THE SPIRAL

Each **process** is a vertical spoke that drives horizontal **MUnit test coverage**.
FeynCalc has 15,222 MUnit tests (~10,000 translatable). Each spiral implements a
new process end-to-end and translates the MUnit tests for the functions it needs.

| Spiral | Process | Status |
|--------|---------|--------|
| 0 | e+e-→μ+μ- tree | DONE (Session 4) |
| 1 | Compton e+γ→e+γ | DONE (Session 5) |
| **2** | **Bhabha e+e-→e+e-** | **NEXT** |
| 3 | QCD qq̄→gg | Planned |
| 4-6 | 1-loop (vertex, VP, Schwinger) | Planned |
| 7 | EW e+e-→W+W- | Planned |
| 8 | MUnit mop-up | Planned |
| 9+ | BSM / ULDM | Planned |

### Branch and code location

- **Branch:** `experimental/rebuild-v2`
- **v2 source:** `src/v2/` (23 files, ~2,700 LOC)
- **v2 tests:** `test/v2/` (9 files, 214 tests)
- **v1:** `src/algebra/`, `src/integrals/` — FROZEN, will be deleted. Do NOT extend or import patterns from.

---

## WHAT EXISTS (v2, 214 tests)

### Six-layer pipeline

```
Layer 1: Model      → qed_model()                    → QEDModel
Layer 2: Rules      → feynman_rules(model)            → FeynmanRules (callable)
Layer 3: Diagrams   → tree_diagrams(model, ...)       → Vector{FeynmanDiagram}
Layer 4: Algebra    → trace → contract → expand → eval → AlgSum (scalar)
Layer 5: Integrals  → PaVe{N}, evaluate(::PaVe; mu2)  → ComplexF64
Layer 6: Evaluate   → solve_tree(prob) → σ             → Float64
```

### Layer 4: Algebra (the core, tree-level complete)

| File | What | Key types/functions |
|------|------|-------------------|
| `coeff.jl` | DimPoly coefficient algebra | `DimPoly`, `DIM`, `Coeff`, `evaluate_dim` |
| `types.jl` | Physics indices, momenta | `LorentzIndex`, `Momentum`, `MomentumSum`, `PairArg` |
| `colour_types.jl` | SU(N) types | `AdjointIndex`, `FundIndex`, `SUNT`, `SUNF`, `SUND`, deltas |
| `pair.jl` | Parametric Pair{A,B} | `MetricTensor`, `FourVector`, `ScalarProduct`, `pair()` |
| `expr.jl` | Dict-based AlgSum | `AlgSum`, `AlgFactor` (6-type union), `FactorKey`, `alg()` |
| `sp_context.jl` | Scalar product context | `SPContext`, `ScopedValues`, `evaluate_sp` |
| `contract.jl` | Lorentz contraction | `contract()`, `substitute_index()`, `alg_from_factors()` |
| `expand_sp.jl` | MomentumSum bilinear expansion | `expand_scalar_product()` |
| `dirac.jl` | Dirac gamma types | `DiracGamma{S}`, `Spinor{K}`, `DiracChain`, `GA/GAD/GS` |
| `dirac_trace.jl` | Trace → AlgSum | `dirac_trace()` (handles up to 8+ gammas) |
| `dirac_expr.jl` | Matrix-valued expressions | `DiracExpr` = Vector{Tuple{AlgSum, DiracChain}} |
| `dirac_trick.jl` | γ^μ...γ_μ contraction | `dirac_trick()` (n=0,1,2 only) |
| `spin_sum.jl` | Fermion spin sums | `spin_sum_amplitude_squared()` (two-line), `_single_line_trace()` |
| `colour_trace.jl` | SU(N) trace | `colour_trace()` (concrete N, recursive) |
| `colour_simplify.jl` | δ contraction, f·f | `contract_colour()` |
| `polarization_sum.jl` | Photon pol sum | `polarization_sum()` (Feynman gauge: -g^{μν}) |

### Layers 1-3, 5-6 (minimal but pipeline-proven)

| File | What |
|------|------|
| `model.jl` | `QEDModel`, `Field{Species}`, `GaugeGroup`, traits |
| `rules.jl` | `FeynmanRules` callable, vertex/propagator dispatch on species |
| `diagrams.jl` | Hard-coded e+e-→μ+μ- s-channel topology |
| `pave.jl` | `PaVe{N}` parametric type, named constructors A0/B0/B1/C0/D0 |
| `pave_eval.jl` | `evaluate(::PaVe{1,2})`, B0 via QuadGK Feynman parameter |
| `schwinger.jl` | Analytical Schwinger correction formula, vacuum polarization |
| `cross_section.jl` | `CrossSectionProblem`, `solve_tree`, Mandelstam, dσ/dΩ |

### Dependencies

- **PolyLog.jl** v2.6.2 — dilogarithm `li2()` (for future C₀ evaluation)
- **QuadGK.jl** v2.11.2 — adaptive quadrature (used by B₀ and vacuum polarization)

### Test files

| File | Tests | What |
|------|-------|------|
| `test_coeff.jl` | 29 | DimPoly arithmetic |
| `test_colour.jl` | 22 | SU(N) traces, δ contraction |
| `test_ee_mumu_x.jl` | 14 | e+e-→μ+μ- algebra (P&S 5.10) |
| `test_self_energy.jl` | 25 | DiracExpr, DiracTrick n=0,1,2 |
| `test_vertical.jl` | 33 | Full pipeline: Model→Rules→Diagrams→Algebra→σ |
| `test_pave.jl` | 51 | PaVe types + A₀/B₀ numerical evaluation |
| `test_schwinger.jl` | 13 | Schwinger correction + vacuum polarization |
| `test_compton.jl` | 4 | Compton |M|² from pipeline vs P&S Eq. 5.87 |
| `test_munit_batch1.jl` | 23 | MUnit translations: DiracTrace, Contract, PolarizationSum |

---

## WHAT WAS DONE IN SESSION 5

### Spiral 0→1 progression

1. **Layer 5 partial:** PaVe{N} types, A₀/B₀ numerical via QuadGK, Schwinger formula
2. **PRD v0.2 rewrite:** Big vision, spiral methodology, anti-hallucination rules, Julia cheatsheet
3. **CLAUDE.md overhaul:** New rules, beads tracking, ground truth acquisition, cheatsheet
4. **Spiral 1 core (Compton):** |M|² matches P&S Eq. 5.87 from algebra pipeline
5. **MUnit batch 1:** 23 translated tests (DiracTrace, Contract, PolarizationSum)
6. **Infrastructure fixes:** AlgSum ==, substitute_index, MomentumSum ordering, expand_sp typing
7. **Beads triage:** 62 closed, 21 open, mapped to spirals

### Key learnings

1. **Don't hand-roll numerics.** First B₀ attempt used analytical formulas → bugs. QuadGK is correct by construction. Rule 8: use Julia packages.

2. **DimD vs Dim4 mismatch is a silent bug.** `GAD(:mu)` creates `LorentzIndex(:mu, DimD())` but `LorentzIndex(:mu)` defaults to `Dim4()`. They print the same but `==` fails. Always use consistent dimension tags.

3. **PaVe is standalone, not in AlgFactor.** PaVe is scalar-valued. It doesn't carry Lorentz indices, doesn't contract, doesn't expand. It does not belong in the AlgFactor union.

4. **substitute_index needed for polarization sums.** When squaring the amplitude, conjugate indices (μ', ν') must be relabeled to (μ, ν) for the contraction engine to find repeated indices. `substitute_index(s, old_li, new_li)` does this.

5. **expand before contract for MomentumSum.** `(p-q)^a (p-q)_a` needs `expand_scalar_product` first (bilinear expansion), then `contract` (repeated index). The other order leaves FourVectors with uncontracted indices.

6. **The FeynCalc Compton example is the ground truth.** File: `refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/ElGa-ElGa.m`. Lines 103-108 give P&S Eq. 5.87 in FeynCalc notation, verified as "CORRECT" by FeynCalc's own test.

---

## WHAT TO DO NEXT: SPIRAL 2 (Bhabha scattering)

### Process: e+e- → e+e- (Bhabha scattering)

This is the next process in the spiral. It requires:

### New capabilities needed

1. **t-channel diagram** — unlike e+e-→μ+μ- (s-channel only) and Compton (s+u), Bhabha has both s-channel AND t-channel. The t-channel has a different topology: the electron exchanges a photon with the positron without annihilating.

2. **DiracTrick n≥3** — Bhabha's squared amplitude involves longer gamma chains. The current DiracTrick only handles n=0,1,2 (γ^μ γ_μ, γ^μ γ^a γ_μ, γ^μ γ^a γ^b γ_μ). Bhabha may need n=3,4. This is the Mertig-Boehm-Denner general formula.
   - MUnit tests: `refs/FeynCalc/Tests/Dirac/DiracTrick.test` — 577 tests, categories `ThreeFreeIndices` (11), `FourFreeIndices` (9), `FiveFreeIndices` (9).

3. **Identical fermion handling** — in Bhabha, the initial and final electrons are the same species. This affects the diagram counting and relative signs (Fermi statistics).

4. **Two-diagram spin sum with interference** — similar to Compton but with s+t channels instead of s+u. The `compton_trace_ij` pattern from test_compton.jl can be reused.

### Ground truth

- **P&S Eq. (5.10) and problem 5.2** — Bhabha in the massless limit
- **FeynCalc example:** look for `refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/ElEl-ElEl.m` or `ElAel-ElAel.m`
- **FeynCalc MUnit:** DiracTrick.test (577 tests), plus Contract.test subset

### Approach

1. Find the FeynCalc Bhabha example in `refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/`
2. Read P&S for the exact formula (Chapter 5, problem 5.2 or nearby)
3. Build the s+t channel amplitudes as gamma chains
4. Compute |M|² using the same trace-contract-evaluate pipeline as Compton
5. Compare to the textbook formula
6. Translate relevant MUnit tests (DiracTrick n≥3 batch)

### Beads issues

Run `bd ready` to see open issues. The relevant ones:
- `feynfeld-7su` — Klein-Nishina (deferred, needs lab-frame kinematics)
- Plus whatever new issues you create for Spiral 2

Create beads issues for Spiral 2 steps before starting implementation:
```bash
bd create --title="Spiral 2.1: Bhabha amplitude (s+t channel)" --description="..." --type=task --priority=1
bd create --title="Spiral 2.2: DiracTrick n>=3" --description="..." --type=task --priority=1
# etc.
```

---

## REFERENCE CODEBASES

All in `refs/` (gitignored):
- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`. Examples in `FeynCalc/Examples/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — Loop integral numerics (Fortran source).
- `refs/papers/` — P&S (.djvu format). More papers needed (Denner 1993, 't Hooft-Veltman) — acquire via TIB VPN.

### Ground truth acquisition

See `CLAUDE.md` §Ground truth acquisition. Priority: arXiv → TIB VPN → playwright-cli.
The P&S textbook is at: `refs/papers/(Frontiers_in_Physics)Michael_E._Peskin,_Dan_V._Schroeder-An_introduction_to_quantum_field_theory-Westview_Press(1995).djvu`
Note: .djvu format. System has `libdjvulibre` but NOT `djvutxt` command-line tool. The `archivum` project at `~/Projects/archivum/` has a djvu extractor but it also needs `djvutxt`. For now, use FeynCalc examples (which cite P&S equations by number) as the citation bridge.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run all v2 tests (214 tests across 9 files)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Run specific test
julia --project=. test/v2/test_compton.jl

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues
bd create --title="..." --description="..." --type=task --priority=1

# Ground truth
ls refs/FeynCalc/Tests/                              # MUnit test directories
ls refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/  # FeynCalc example computations
ls refs/papers/                                      # local paper copies

# Commit and push (session end protocol)
git add <files>
git commit -m "..."
git push
bd backup export-git  # sync beads to origin
```

---

## FILE MAP

```
src/v2/
├── FeynfeldX.jl          # Module root, includes + exports (89 LOC)
├── coeff.jl              # DimPoly coefficients (142)
├── types.jl              # PhysicsIndex, Momentum, MomentumSum (94)
├── colour_types.jl       # SU(N): AdjointIndex, FundIndex, SUNT, SUNF, SUND (126)
├── pair.jl               # Parametric Pair{A,B}, NOT exported (76)
├── expr.jl               # AlgSum (Dict), AlgFactor (6-type union), FactorKey (151)
├── sp_context.jl         # SPContext + ScopedValues (71)
├── contract.jl           # Lorentz contraction + substitute_index (155)
├── expand_sp.jl          # Scalar product bilinear expansion (83)
├── dirac.jl              # DiracGamma{S}, Spinor{K}, DiracChain (113)
├── dirac_trace.jl        # Dirac trace → AlgSum (65)
├── dirac_expr.jl         # DiracExpr: matrix-valued expressions (104)
├── dirac_trick.jl        # D-dim γ^μ...γ_μ for n=0,1,2 (117)
├── spin_sum.jl           # Fermion spin sums (139)
├── colour_trace.jl       # SU(N) trace → AlgSum (82)
├── colour_simplify.jl    # Delta contraction, f·f identity (148)
├── polarization_sum.jl   # Photon pol sum: -g^{μν} (16)
├── model.jl              # AbstractModel, QEDModel, Field{Species} (101)
├── rules.jl              # FeynmanRules callable (81)
├── diagrams.jl           # FeynmanDiagram, hard-coded topologies (75)
├── pave.jl               # PaVe{N} type, named constructors (68)
├── pave_eval.jl          # evaluate(::PaVe) via QuadGK (90)
├── schwinger.jl          # Schwinger correction + vacuum pol (53)
├── cross_section.jl      # Mandelstam, Problem/Solve, σ (108)
├── DESIGN.md             # Design choices, anti-patterns, cockroaches
└── VERTICAL_PLAN.md      # Original vertical plan (historical)

test/v2/
├── test_coeff.jl         # DimPoly (29 tests)
├── test_colour.jl        # SU(N) (22 tests)
├── test_ee_mumu_x.jl     # e+e-→μ+μ- algebra (14 tests)
├── test_self_energy.jl   # DiracExpr + DiracTrick (25 tests)
├── test_vertical.jl      # Full pipeline (33 tests)
├── test_pave.jl          # PaVe types + numerics (51 tests)
├── test_schwinger.jl     # Schwinger correction (13 tests)
├── test_compton.jl       # Compton |M|² vs P&S 5.87 (4 tests)
└── test_munit_batch1.jl  # MUnit translations (23 tests)

Total: ~2,700 source LOC, ~900 test LOC, 214 tests, all files < 200 LOC.
```
