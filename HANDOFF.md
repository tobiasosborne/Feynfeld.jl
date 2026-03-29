# HANDOFF — 2026-03-29 (End of Session 9, Spiral 9 + P1 bugs complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Read `Feynfeld_PRD.md` — vision, endgame interaction, spiral plan, hardness scale
3. Read `src/v2/DESIGN.md` — type system, anti-patterns, Session 8 review findings
4. Run `bd ready` to see available work
5. Run `julia --project=. test/v2/test_vertical.jl` to verify pipeline works
6. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## THE PIPELINE PRINCIPLE

**The pipeline IS the architecture. If a process bypasses it, the architecture is incomplete.**

Every process MUST flow through: Model → Rules → Channels → Amplitude → Algebra → Evaluate.
No hand-built amplitudes in test files. No standalone recipes that skip the pipeline.
Standalone formulas (schwinger.jl, ew_cross_section.jl, etc.) are REFERENCE IMPLEMENTATIONS
for cross-validation — not substitutes for the pipeline.

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **CITE EVERYTHING (tiered).**
3. **ACQUIRE PAPERS BEFORE CODE.** arXiv → TIB VPN → playwright-cli.
4. **JULIA IDIOMATIC.** Dispatch, not isa cascades. No Any.
5. **WORKFLOW (TIERED):** <5 LOC direct, <20 LOC 1+1, >20 LOC 3+1.
6. **REVIEW.** Rigorous reviewer after every core change.
7. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.
8. **LOC LIMIT ~200.** No source file exceeds ~200 lines.
9. **NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT OVERVIEW

### What is Feynfeld.jl?

Julia-native, agent-facing, full-stack physics computation suite. Lagrangian →
cross-section in one `using Feynfeld`. Replaces FeynRules + FeynArts + FeynCalc +
FormCalc + LoopTools. See PRD for the full vision.

### Branch and code location

- **Branch:** `experimental/rebuild-v2`
- **v2 source:** `src/v2/` (34 files, ~4,200 LOC)
- **v2 tests:** `test/v2/` (17 files, 329 tests)
- **v1:** `src/algebra/`, `src/integrals/` — FROZEN, will be deleted. Do NOT extend.

---

## WHAT EXISTS (v2, 329 tests)

### Six-layer pipeline (now with channel enumeration)

```
Layer 1: Model      → qed_model() / qcd_model() / ew_model()  → AbstractModel
Layer 2: Rules      → feynman_rules(model)                      → FeynmanRules
Layer 3: Channels   → tree_channels(model, rules, in, out)      → Vector{TreeChannel}
         Amplitude  → build_amplitude(ch, rules, model)          → DiracChains
Layer 4: Algebra    → trace → contract → expand → eval           → AlgSum (scalar)
Layer 5: Integrals  → PaVe{N}, evaluate(::PaVe; mu2)            → ComplexF64
Layer 6: Evaluate   → solve_tree(prob) → σ                      → Float64
```

### Pipeline coverage (Session 9 result)

| Process | Channels | Status |
|---------|----------|--------|
| e+e-→μ+μ- | s(γ) | Full pipeline via solve_tree |
| Bhabha | s(γ)+t(γ) | Pipeline + spin_sum_interference |
| Compton | s(e)+u(e) | Pipeline + _cross_line_trace |
| qq̄→gg | t(q)+u(q)+s(g) | Pipeline t/u, manual s (ggg vertex) |
| ee→WW | s(γ)+s(Z)+t(ν)+u(ν) | Channel enumeration only |

### Source files

| File | LOC | What |
|------|-----|------|
| **Layer 4: Algebra** | | |
| `coeff.jl` | 142 | DimPoly coefficient algebra |
| `types.jl` | 79 | LorentzIndex, Momentum, MomentumSum |
| `colour_types.jl` | 126 | SUNT, SUNF, SUND, deltas |
| `pair.jl` | 76 | Parametric Pair{A,B} |
| `expr.jl` | 149 | AlgSum (Dict), AlgFactor, FactorKey |
| `sp_context.jl` | 71 | SPContext + ScopedValues |
| `contract.jl` | 136 | Lorentz contraction + Eps handling |
| `eps_contract.jl` | 85 | ε·ε = -det[pair(aᵢ,bⱼ)] |
| `expand_sp.jl` | 82 | Scalar product bilinear expansion |
| `dirac.jl` | 118 | DiracGamma{S}, Spinor{K}, DiracChain |
| `dirac_trace.jl` | ~160 | Trace → AlgSum (gamma5, projectors) |
| `dirac_expr.jl` | 104 | DiracExpr: matrix-valued expressions |
| `dirac_trick.jl` | 117 | D-dim γ^μ...γ_μ for n=0..5+ |
| `spin_sum.jl` | 138 | Fermion spin sums (completeness) |
| `interference.jl` | 98 | Cross-line traces, spin_sum_interference |
| `colour_trace.jl` | 72 | SU(N) trace → (real, imag) AlgSum (**FIXED**: i² on f-terms) |
| `colour_simplify.jl` | 148 | Delta contraction via dispatch |
| `polarization_sum.jl` | 58 | Feynman/axial/massive pol sums |
| **Layers 1-3** | | |
| `model.jl` | 99 | AbstractModel, QEDModel, Field{Species} |
| `qcd_model.jl` | 60 | QCDModel, qqg + ggg vertices, triple_gauge_vertex |
| `ew_model.jl` | 49 | EWModel, 5 SM vertex types |
| `rules.jl` | 82 | FeynmanRules callable, vertex dispatch |
| `diagrams.jl` | 22 | ExternalLeg (mass field, backward-compat) |
| `channels.jl` | 103 | TreeChannel, tree_channels() |
| `amplitude.jl` | 141 | build_amplitude: boson + fermion exchange |
| **Layer 5: Integrals** | | |
| `pave.jl` | ~80 | PaVe{N} type, named constructors |
| `pave_eval.jl` | 182 | evaluate: A₀ closed, B₀ hybrid (quadgk real + Kallen imag), C₀/C₁/C₂ quadgk (**SLOW**) |
| **Layer 6 + Reference** | | |
| `cross_section.jl` | ~100 | Mandelstam, solve_tree, σ |
| `schwinger.jl` | ~50 | REFERENCE: Schwinger correction |
| `vertex.jl` | 69 | REFERENCE: QED g-2 F₂(0)=α/(2π) |
| `running_alpha.jl` | ~100 | REFERENCE: running α(q²) |
| `ew_parameters.jl` | 33 | EW constants: M_W, M_Z, sin²θ_W |
| `ew_cross_section.jl` | 83 | REFERENCE: σ(ee→WW) Grozin formula |

### Test files (329 tests across 17 files)

| File | Tests | What |
|------|-------|------|
| `test_coeff.jl` | 29 | DimPoly arithmetic |
| `test_colour.jl` | 27 | SU(N) traces (**incl n≥4 fixed**), δ contraction, f·f/d·d |
| `test_ee_mumu_x.jl` | 14 | e+e-→μ+μ- algebra (P&S 5.10) |
| `test_self_energy.jl` | 25 | DiracExpr, DiracTrick n=0,1,2 |
| `test_vertical.jl` | 34 | Full pipeline via solve_tree |
| `test_pave.jl` | 53 | PaVe types + A₀/B₀/B₁ (**incl Im(B₀) analytical**) |
| `test_schwinger.jl` | 15 | Schwinger correction + vacuum polarization |
| `test_compton.jl` | 4 | Compton |M|² from pipeline vs P&S 5.87 |
| `test_munit_batch1.jl` | 23 | MUnit: DiracTrace, Contract, PolarizationSum |
| `test_munit_batch2.jl` | 18 | MUnit: DiracTrick n=3,4 |
| `test_bhabha.jl` | 4 | Bhabha |M̄|² at 2 kinematic points |
| `test_qqbar_gg.jl` | 2 | QCD qq̄→gg |M̄|² at 2 kinematic points |
| `test_self_energy_1loop.jl` | 13 | 1-loop Σ(p) (**assertions fixed for correct Im(B₀)**) |
| `test_vertex_g2.jl` | 32 | C₀/C₁/C₂, F₂(0)=α/(2π) (**SLOW: 8min, needs analytical C₀**) |
| `test_running_alpha.jl` | 34 | Running α(q²), Δα, improved Born σ |
| `test_ee_ww.jl` | 36 | Tree-level e⁺e⁻→W⁺W⁻ reference formula |
| `test_pipeline.jl` | 17 | **NEW**: Bhabha/Compton/qq→gg/ee→WW pipeline tests |

---

## WHAT WAS DONE IN SESSION 9

### Spiral 9: Pipeline consolidation (Phases A-D)

**6 beads issues closed.** New files: channels.jl, amplitude.jl, interference.jl,
qcd_model.jl, ew_model.jl, test_pipeline.jl.

- **Phase A**: TreeChannel + tree_channels() + build_amplitude for boson exchange
- **Phase B**: Fermion exchange + spin_sum_interference (Bhabha + Compton)
- **Phase C**: QCDModel with qqg/ggg vertices, qq→gg pipeline test
- **Phase D**: EWModel with 5 SM vertices, channel enumeration for ee→WW

### P1 bug fixes

- **feynfeld-83m (colour_trace i²)**: Recursive trace now returns (real, imag)
  tuple. The product i·f × i·f = -f·f is correctly handled. Verified:
  Tr(T^aT^bT^bT^a) = 16/3, Tr(T^aT^bT^aT^b) = -2/3 for N=3.

- **feynfeld-1rb (B₀ imaginary part)**: Analytical Im(B₀) = π√λ/p² via
  Kallen function λ(p²,m₀²,m₁²). One-massless case: Im = π(p²-m²)/p².
  test_self_energy_1loop assertions updated (were masked by the bug).

### Cleanup

- spin_sum.jl split → spin_sum.jl (138) + interference.jl (98)
- 4 standalone recipe files marked as REFERENCE IMPLEMENTATION
- ExternalLeg gains mass field (backward-compatible default 0//1)
- All amplitude indices use DimD() for D-dimensional traces
- evaluate_m_squared calls evaluate_dim for DimPoly coefficients

---

## WHAT TO DO NEXT

### Priority 1: Performance — analytical loop integrals

**THE #1 BOTTLENECK.** test_vertex_g2 takes 8 minutes because C₀ uses nested
quadgk (O(N²) function evaluations). The field solved this in 1979.

| Issue | What | Impact |
|-------|------|--------|
| `feynfeld-fqy` | Analytical C₀ via 't Hooft-Veltman Spence formulas | **8min → <1sec** |
| `feynfeld-0ku` | Analytical B₀ closed-form (logs only) | ~2x speedup |
| `feynfeld-2xv` | Analytical D₀ (needed for box diagrams, Spiral 10) | Unblocks Spiral 10 |

**How to implement C₀:**
1. Read `refs/papers/tHooftVeltman1979_NuclPhysB153.pdf` Eqs (5.2)-(5.30)
2. Read `refs/LoopTools/src/C/C0func.F` for the case dispatch pattern
3. Julia already has `PolyLog.jl` with `li2` (dilogarithm) — use it
4. Dispatch on number of zero momenta: C0p0/C0p1/C0p2/C0p3
5. Handle soft/collinear limits as separate cases
6. ~200-400 LOC, all in a new `pave_analytical.jl`
7. Keep quadgk as numerical fallback for validation

**Expected result:**
```
C₀ evaluation: 1-5 ms → 1-10 μs (1000x speedup)
test_vertex_g2: 8 min → < 1 sec
Full test suite: 12 min → < 2 min
```

### Priority 2: Test runner performance

| Issue | What | Impact |
|-------|------|--------|
| `feynfeld-icg` | Single-process test runner (1 Julia process, not 17) | **60s → 5s** JIT |
| `feynfeld-6nu` | Make FeynfeldX precompilable (blocked by rename) | Near-instant startup |

**How to implement test runner:**
Create `test/v2/runtests.jl`:
```julia
using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX
@testset "FeynfeldX" begin
    include("test_coeff.jl")  # each file uses @testset but not include/using
    include("test_colour.jl")
    ...
end
```
Each test file needs the `include/using` lines removed (guarded by `@isdefined`).

### Priority 3: Remaining pipeline work

| Issue | What |
|-------|------|
| ee→WW full pipeline | Chiral vertices (eeZ with γ5), triple gauge (WWγ/WWZ), massive propagators |
| `feynfeld-qyu` | Rename FeynfeldX → Feynfeld |
| `feynfeld-akr` | CI/CD GitHub Actions |

### Priority 4: Spiral 10 — D₀ + box diagrams

Requires analytical D₀ (feynfeld-2xv) which requires analytical C₀ (feynfeld-fqy).
The dependency chain: C₀ analytical → D₀ analytical → box diagrams → NLO physics.

---

## KNOWN BUGS (all P1 fixed, remaining are P2+)

- **vertex_f2 above threshold** — doesn't implement iε for q² > 4m². Use only
  for spacelike or below threshold.
- **B₀ real part** still uses quadgk (correct but slow). feynfeld-0ku tracks this.
- **C₀/D₀** use nested quadgk (correct but 1000x slow). feynfeld-fqy/2xv track this.

---

## REFERENCE CODEBASES

All in `refs/` (gitignored):
- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — **KEY for analytical integrals**: C₀ in `src/C/`, D₀ in `src/D/`.
- `refs/papers/` — 16+ local paper copies.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run specific test (fast)
julia --project=. test/v2/test_vertical.jl    # 5s, pipeline
julia --project=. test/v2/test_pipeline.jl     # 7s, all processes
julia --project=. test/v2/test_colour.jl       # 3s, incl n>=4

# Fast smoke test (~40s, skips slow vertex_g2)
for f in test/v2/test_coeff.jl test/v2/test_colour.jl test/v2/test_vertical.jl \
         test/v2/test_pipeline.jl test/v2/test_pave.jl; do
    julia --project=. "$f"
done

# Full suite (~12 min, vertex_g2 is 8min alone — fix with analytical C₀)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```

---

## FILE MAP

```
src/v2/
├── FeynfeldX.jl          # Module root, includes + exports
├── coeff.jl              # DimPoly coefficients
├── types.jl              # PhysicsIndex, Momentum, MomentumSum
├── colour_types.jl       # SU(N): AdjointIndex, FundIndex, SUNT, SUNF, SUND
├── pair.jl               # Parametric Pair{A,B}
├── expr.jl               # AlgSum (Dict), AlgFactor, FactorKey
├── sp_context.jl         # SPContext + ScopedValues
├── contract.jl           # Lorentz contraction + Eps index handling
├── eps_contract.jl       # ε·ε → -det[pair(aᵢ,bⱼ)]
├── expand_sp.jl          # Scalar product bilinear expansion
├── dirac.jl              # DiracGamma{S}, Spinor{K}, DiracChain
├── dirac_trace.jl        # Dirac trace (gamma5 + projectors)
├── spin_sum.jl           # Fermion spin sums (completeness)
├── interference.jl       # NEW: cross-line traces, reconnected interference
├── dirac_expr.jl         # DiracExpr: matrix-valued expressions
├── dirac_trick.jl        # D-dim γ^μ...γ_μ
├── colour_trace.jl       # FIXED: SU(N) trace with (real, imag) tracking
├── colour_simplify.jl    # Delta contraction via dispatch
├── polarization_sum.jl   # Feynman + axial + massive pol sums
├── model.jl              # AbstractModel, QEDModel
├── qcd_model.jl          # NEW: QCDModel, triple_gauge_vertex
├── ew_model.jl           # NEW: EWModel with 5 SM vertices
├── rules.jl              # FeynmanRules callable
├── diagrams.jl           # ExternalLeg (with mass field)
├── channels.jl           # NEW: TreeChannel, tree_channels()
├── amplitude.jl          # NEW: build_amplitude (boson + fermion exchange)
├── pave.jl               # PaVe{N} type, named constructors
├── pave_eval.jl          # FIXED: B₀ imag via Kallen. SLOW: C₀ still quadgk
├── schwinger.jl          # REFERENCE: Schwinger correction
├── vertex.jl             # REFERENCE: QED g-2 F₂(0)
├── running_alpha.jl      # REFERENCE: SM running α(q²)
├── ew_parameters.jl      # EW: M_W, M_Z, sin²θ_W, Z couplings
├── ew_cross_section.jl   # REFERENCE: σ(ee→WW) Grozin formula
├── cross_section.jl      # Mandelstam, solve_tree, σ
└── DESIGN.md             # Design choices + review findings

test/v2/                  # 17 files, 329 tests
reviews/                  # 6 review reports from Session 8
SPIRAL_9_PLAN.md          # Spiral 9 plan (completed)
```
