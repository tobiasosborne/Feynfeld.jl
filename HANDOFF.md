# HANDOFF — 2026-03-30 (Session 12: D₀ fixed, EW pipeline, Spiral 10 started)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (396+ tests, ~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 12 ACCOMPLISHMENTS

### 1. Fixed the D₀ compilation bomb
**Root cause:** Session 11 agent's `_D0_quadgk` had triple-nested `quadgk` closures. Julia JIT
compiles ALL reachable code paths, including unreachable fallbacks. Triple-nested closures
cause a compilation time explosion (minutes, gigabytes of RAM).

**Fix:** Removed `_D0_quadgk` fallback. `_D0_evaluate` now errors if COLLIER unavailable.
COLLIER works correctly — returns in microseconds.

### 2. Removed LoopTools dependency
LoopTools was declared in Project.toml but NEVER imported in any source file. Its precompilation
(Fortran FFI wrapping) added minutes to startup. Removed from `[deps]` and `[compat]`.

### 3. Single-process test runner
`test/v2/runtests.jl` — loads FeynfeldX once, includes all 18 test files. Each test file
has `@isdefined(FeynfeldX) || include(...)` guard for standalone use. 396 tests in 5m15s.

### 4. Pipeline completion for all 5 processes
| Process | Pipeline Status |
|---------|----------------|
| e+e-→μ+μ- | Full ✓ |
| Bhabha | Full ✓ |
| Compton | Full ✓ |
| qq̄→gg | **Full ✓** (gauge exchange dispatch added) |
| ee→WW | **Builds all channels ✓**, Grozin comparison pending |

### 5. Vertex structure dispatch (EW infrastructure)
`vertex_structure` now dispatches on `Val(coupling)`:
- `Val(:e)` → γ^μ (QED/QCD)
- `Val(:e_Z)` → (g_V - g_A γ5) γ^μ (neutral current, uses Rational EW_GV_E_R)
- `Val(:g_W)` → (1-γ5)/2 γ^μ (charged current)
- `Val(:g_s)` → γ^μ (QCD)

`build_amplitude` now uses `_lookup_vertex` to get vertex structures from FeynmanRules.
Returns `DiracExpr` (not DiracChain) to support chiral vertices with multiple terms.
`spin_sum_amplitude_squared` extended for DiracExpr. Dirac conjugation: GA6↔GA7 swap.

### 6. D-tensor PV reduction (Spiral 10)
`d_tensor.jl` (57 LOC): D₁/D₂/D₃ via Passarino-Veltman reduction with 3×3 Gram matrix.
4 C₀ sub-integrals (one per removed propagator). Symmetric check: D₁=D₂=D₃ to 1e-12.

### 7. MUnit test porting started
`test/v2/munit/test_DiracTrace.jl`: 22 symbolic tests from FeynCalc DiracTrace.test.
All comparisons are exact AlgSum equality — NO numerical spot-checks.

---

## KNOWN ISSUES AND BLOCKERS

### P1: Rational{Int} overflow in AlgSum (feynfeld-6ds)
The coefficient system uses `Rational{Int}` which overflows when Float64 kinematics are
rationalized and multiplied with EW coupling constants. Blocks the full ee→WW Grozin
cross-section comparison. Fix options:
1. BigRational fallback in coeff.jl
2. Float64 coefficient path
3. Pre-contract symbolic expressions to Float64 evaluation functions

### ee→WW Grozin comparison (feynfeld-bao)
All 3 physical channels (s-γ, s-Z, t-ν) build correctly through the pipeline.
Diagonal |M|² evaluates at integer kinematic points. Full cross-section integration
blocked by Rational overflow during quadgk (non-integer cosθ values).

---

## WHAT TO DO NEXT

### Priority 1: MUnit test porting (~386 remaining tests)

**12 beads created** covering ~408 portable tests. 22 done, 386 remaining.

**Immediate next:** Finish DiracTrace (36 tests remaining). The big ones:
- ID41-42: 8-gamma traces (105 terms each) — **need programmatic Mathematica→Julia translator**
- ID43-44: 6-gamma + γ5 → ε·SP mixed terms
- ID56-59: 6/8-gamma + γ5 → ε·MTD mixed terms (15-105 terms)
- ID14-16, ID20, ID24, ID28: projector + momentum combos

**CRITICAL RULE:** All test comparisons MUST be exact symbolic (AlgSum ==).
Numerical spot-checks are "corruption and poison." Even 105-term expressions must be
translated programmatically into full symbolic expected values.

**Translator approach:** Write a script that converts Mathematica `SP[a,b]*SP[c,d]` → Julia
`alg(SP(:a,:b))*alg(SP(:c,:d))`. The notation is mechanical:
- `SP[a, b]` → `alg(SP(:a, :b))`
- `MT[i, j]` → `alg(MT(:i, :j))` (4D) or `alg(MTD(:i, :j))` (D-dim)
- `LCD[a,b,c,d]` → `alg(Eps(LorentzIndex(:a,DimD()), ...))`
- `LC[a,b,c,d]` → `alg(Eps(LorentzIndex(:a), ...))`
- Multiplication `*` stays `*`, addition `+`/`-` stays

**Skippable IDs** (not implementable):
- ID1-6: unevaluated structural (not computational)
- ID29-31, ID46: scheme-dependent (BMHV/NDR)
- ID33, ID53: DiracSigma (not implemented)
- ID34-39: Cartesian gammas (CGA/CGS — not supported)
- ID47-51: DOT/scalar multiplication (not trace tests)

**Dimension convention:** `GAD(:i)` creates DimD() indices. Expected values must use
`MTD(:i,:j)` (not `MT(:i,:j)` which creates Dim4). `SP(:p,:q)` is dimension-agnostic.

**MUnit bead IDs:**
| Issue | Category | Tests |
|-------|----------|-------|
| feynfeld-q6m | DiracTrace | 58 (22 done) |
| feynfeld-32j | DiracTrick 2-idx | 27 |
| feynfeld-37v | DiracTrick 3-idx | 11 |
| feynfeld-iaz | DiracTrick 4-idx | 9 |
| feynfeld-n01 | DiracTrick 5-idx | 9 |
| feynfeld-8qe | DiracTrick 1-idx | 17 |
| feynfeld-s4p | EpsContract | 41 |
| feynfeld-rcd | Contract | 20 |
| feynfeld-8xb | ExpandScalarProduct | 15 |
| feynfeld-36h | SUNTrace | 24 |
| feynfeld-4mm | SUNSimplify | 78 |
| feynfeld-mfb | PolarizationSum | ~15 |

### Priority 2: Spiral 10 continuation

**D₀/D₁/D₂/D₃ evaluation: DONE.** All pass tests.

**Next steps:**
- `feynfeld-7h8`: 1-loop amplitude builder for 2→2 box diagrams (~100-150 LOC)
- `feynfeld-4q5`: ee→μμ NLO box diagram via pipeline (blocked by builder)

The box diagram for ee→μμ has 1 QED photon loop. 4 propagators forming a box.
After PaVe decomposition: D₀ + tensor coefficients. Validate against known NLO correction.

### Priority 3: Rational overflow fix (feynfeld-6ds)
Unblocks ee→WW Grozin comparison. Affects any computation mixing Float64 kinematics
with Rational symbolic coefficients.

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
8. **NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT STATE

### Branch and code location
- **Branch:** `master`
- **v2 source:** `src/v2/` (37 files, ~3,900 LOC)
- **v2 tests:** `test/v2/` (19 files + munit/) — 396+ tests
- **v1:** FROZEN. Do not extend.

### Source files (37 files)

| File | LOC | What |
|------|-----|------|
| **Layer 4: Algebra** | | |
| `coeff.jl` | 142 | DimPoly coefficient algebra |
| `types.jl` | 79 | LorentzIndex, Momentum, MomentumSum |
| `colour_types.jl` | 126 | SUNT, SUNF, SUND, deltas |
| `pair.jl` | 76 | Parametric Pair{A,B}, SP/MT/MTD helpers |
| `expr.jl` | 149 | AlgSum (Dict), AlgFactor, FactorKey |
| `sp_context.jl` | 71 | SPContext + ScopedValues |
| `contract.jl` | 136 | Lorentz contraction + Eps handling |
| `eps_contract.jl` | 85 | ε·ε = -det[pair(aᵢ,bⱼ)] |
| `expand_sp.jl` | 82 | Scalar product bilinear expansion |
| `dirac.jl` | 118 | DiracGamma{S}, Spinor{K}, DiracChain |
| `dirac_trace.jl` | ~160 | Trace → AlgSum (gamma5, projectors) |
| `dirac_expr.jl` | 104 | DiracExpr: matrix-valued expressions |
| `dirac_trick.jl` | 117 | D-dim γ^μ...γ_μ for n=0..5+ |
| `spin_sum.jl` | ~170 | Fermion spin sums (DiracChain + DiracExpr) |
| `interference.jl` | ~105 | Cross-line traces, spin_sum_interference |
| `colour_trace.jl` | 72 | SU(N) trace → (real, imag) AlgSum |
| `colour_simplify.jl` | 148 | Delta contraction via dispatch |
| `polarization_sum.jl` | 58 | Feynman/axial/massive pol sums |
| **Layers 1-3** | | |
| `model.jl` | 99 | AbstractModel, QEDModel, Field{Species} |
| `qcd_model.jl` | ~70 | QCDModel, qqg + ggg vertices, _MomLike |
| `ew_model.jl` | 49 | EWModel, 5 SM vertex types |
| `rules.jl` | ~95 | FeynmanRules, vertex_structure w/ Val(coupling) |
| `diagrams.jl` | 22 | ExternalLeg (mass field) |
| `channels.jl` | 103 | TreeChannel, tree_channels() |
| `amplitude.jl` | ~180 | build_amplitude: boson/fermion/gauge exchange |
| **Layer 5: Integrals** | | |
| `pave.jl` | ~80 | PaVe{N} type, A0/B0/C0/D0/D1/D2/D3 constructors |
| `pave_eval.jl` | ~215 | evaluate: A₀/B₀/B₁ + C₁/C₂ PV reduction |
| `c0_analytical.jl` | 82 | C₀: COLLIER ccall + C0p0 analytical |
| `d0_collier.jl` | 42 | D₀: COLLIER ccall (no quadgk fallback) |
| `d_tensor.jl` | 57 | D₁/D₂/D₃: PV reduction, 3×3 Gram matrix |
| **Layer 6 + Reference** | | |
| `cross_section.jl` | ~100 | Mandelstam, solve_tree, σ |
| `schwinger.jl` | ~50 | REFERENCE: Schwinger correction |
| `vertex.jl` | 69 | REFERENCE: QED g-2 F₂(0)=α/(2π) |
| `running_alpha.jl` | ~100 | REFERENCE: running α(q²) |
| `ew_parameters.jl` | ~40 | EW constants (Float64 + Rational versions) |
| `ew_cross_section.jl` | 83 | REFERENCE: σ(ee→WW) Grozin formula |

---

## COLLIER SETUP (REQUIRED ON EACH MACHINE)

```bash
cd refs/COLLIER/COLLIER-1.2.8
mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd /path/to/Feynfeld.jl
ls refs/COLLIER/COLLIER-1.2.8/libcollier.so  # verify
mkdir -p output  # COLLIER writes log files here
```

---

## QUICK COMMANDS

```bash
# Run single test (fast)
julia --project=. test/v2/test_vertical.jl    # 5s, pipeline
julia --project=. test/v2/test_d0.jl          # 30s, D₀+D-tensor
julia --project=. test/v2/munit/test_DiracTrace.jl  # 2s, MUnit

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
