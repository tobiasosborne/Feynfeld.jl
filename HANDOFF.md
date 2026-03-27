# HANDOFF — 2026-03-27 (End of Session 4)

## DO NOT DELETE THIS FILE. Read it completely before working.

## TOBIAS'S RULES — FOLLOW TO THE LETTER

1. **SKEPTICISM**: All subagent work, handoffs — verify everything twice.
2. **DEEP BUGS**: Deep, complex, interlocked. Do not underestimate.
3. **NO BANDAIDS**: Best-practices full solutions only.
4. **WORKFLOW**: 3 subagents before any core code change (research source + 2 solutions).
5. **REVIEW**: Rigorous reviewer agent after every core change. No exceptions.
6. **GROUND TRUTH**: Physics is ground truth, not pinned numbers. Tests may be suspect.
7. **TESTING**: Targeted only, or full suite in background.
8. **REPEAT RULES**: Repeat occasionally to maintain focus.
9. **DO NOT UNDERESTIMATE**: This is deeply nontrivial.
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Read-only research/design agents CAN run in parallel.
11. **RESEARCH JULIA IDIOMS FIRST**: Before every new layer, research Julia-idiomatic patterns. Don't default to OOP struct hierarchies. (Added Session 4.)

**NEVER modify TensorGR.jl without explicit permission.** It is an active separate
project with its own workflow and handoff protocol.

---

## CRITICAL CONTEXT: v1 vs v2

**There are TWO implementations. v2 is the future.**

### v1 (master branch, `src/algebra/`, `src/integrals/`)
- FeynCalc-mirror design: tagged unions, `Expr` coefficients, polymorphic returns
- Phase 1 algebra COMPLETE: 616 tests pass
- Has: PaVe types, FAD, Tdec, PaVeReduce (integrals layer), TensorGR bridge
- Does NOT have: Model, Rules, Diagrams, cross-section evaluation
- **Status: FROZEN. Do not extend. Use only as reference for porting algorithms.**

### v2 (branch `experimental/rebuild-v2`, `src/v2/`)
- Julia-idiomatic: parametric types, DimPoly coefficients, Dict-based AlgSum
- Full pipeline validated: Model → Rules → Diagrams → Algebra → Evaluate
- 123 tests pass, 2,056 source LOC across 20 files
- **Status: ACTIVE. All new work goes here.**

**Read `src/v2/DESIGN.md` before touching any v2 code.** It documents every design
choice, every cockroach found, and the patterns to follow.

---

## What's Done (v2)

### Algebra (Layer 4) — COMPLETE for tree-level + 1-loop numerator

| Component | File | Tests | Status |
|-----------|------|-------|--------|
| DimPoly coefficients | `coeff.jl` | 29 | D, D-4, D², arithmetic, normalisation |
| Core types | `types.jl` | — | PhysicsIndex, LorentzIndex, Momentum, MomentumSum |
| Parametric Pair{A,B} | `pair.jl` | — | MetricTensor/FourVector/ScalarProduct aliases |
| Dict-based AlgSum | `expr.jl` | — | 6-type AlgFactor union, FactorKey, full arithmetic |
| SPContext + ScopedValues | `sp_context.jl` | — | Implicit via ScopedValues, explicit override |
| Lorentz contraction | `contract.jl` | — | Worklist with dispatch on Pair types |
| SP expansion | `expand_sp.jl` | — | Bilinear MomentumSum expansion |
| DiracGamma{S} | `dirac.jl` | — | Parametric slots, Spinor{K}, DiracChain |
| Dirac trace | `dirac_trace.jl` | 14 | Always returns AlgSum, MomSumSlot expansion |
| DiracExpr | `dirac_expr.jl` | 25 | Matrix-valued expressions, arithmetic, simplify |
| DiracTrick | `dirac_trick.jl` | 25 | D-dim γ^μ...γ_μ for n=0,1,2 |
| Spin sum | `spin_sum.jl` | 14 | Separate traces per fermion line, multiply |
| SU(N) colour types | `colour_types.jl` | — | AdjointIndex, FundIndex, SUNT, SUNF, SUND, deltas |
| Colour trace | `colour_trace.jl` | 22 | Concrete N, recursive for n≥4 |
| Colour simplify | `colour_simplify.jl` | 22 | δ contraction, f·f = Nδ identity |

### Higher Layers — MINIMAL but pipeline-proven

| Layer | File | Tests | Status |
|-------|------|-------|--------|
| 1: Model | `model.jl` | 13 | AbstractModel interface, Field{Species}, GaugeGroup types, traits |
| 2: Rules | `rules.jl` | 7 | Callable FeynmanRules, dispatch on species for propagators/vertices |
| 3: Diagrams | `diagrams.jl` | 6 | Hard-coded e+e-→μ+μ- s-channel; full generation NOT implemented |
| 6: Evaluate | `cross_section.jl` | 3 | Mandelstam, Problem pattern, σ_total analytical |

### Validation

| Test | What | Result |
|------|------|--------|
| `test_ee_mumu_x.jl` | Algebra tracer: P&S Eq (5.10) | |M|² = 8(t²+u²) ✓ |
| `test_coeff.jl` | DimPoly: D², D*D, normalisation | All 29 ✓ |
| `test_colour.jl` | SU(N): Tr(T^aT^b), δ^{aa}, f·f contraction | All 22 ✓ |
| `test_self_energy.jl` | DiracTrick: γ^μ(m+γ·q)γ_μ = Dm+(2-D)γ·q | All 25 ✓ |
| `test_vertical.jl` | **FULL PIPELINE**: Model→Rules→Diagrams→Algebra→σ | P&S (5.10) + (5.12) ✓ |

---

## What's NOT Done

### Layer 5: Integrals (the big gap)

v2 has NO loop integral machinery. v1 has symbolic PaVe types that can be ported.

Needed for 1-loop calculations:
1. **PaVe types** — port `PaVe{N}` from v1, adapt to v2 DimPoly coefficients
2. **Tensor decomposition** — port `tdec()` from v1 for rank 0-2
3. **PaVe reduction** — port B-function reduction from v1
4. **Numerical evaluation** — NEW: Li₂, A₀, B₀ in pure Julia (no Fortran FFI)
5. **ToPaVe** — FAD → PaVe extraction (v1 doesn't have this either)

### Layer 3: Diagrams (hard-coded only)

Only e+e- → μ+μ- s-channel is hard-coded. Full topology generation not implemented.
For the next process (Compton, Bhabha, etc.), either hard-code or implement:
- Recursive tree topology enumeration
- Field insertion at vertices
- Momentum routing

### Gaps in existing layers

| Gap | Layer | Impact | Difficulty |
|-----|-------|--------|-----------|
| DiracTrick n≥3 | Algebra | 1-loop with ≥3 internal propagators | Medium — Mertig/Boehm/Denner formula |
| γ5 traces n>4 | Algebra | Axial-vector processes | Medium |
| BMHV/Larin scheme | Algebra | 2-loop γ5 | Hard |
| Fierz identity | Algebra | Multi-loop colour | Medium |
| Tdec rank ≥3 | Algebra | Box diagrams | Medium |
| PaVeReduce C/D | Integrals | 3/4-point functions | Hard — Denner recursion |
| Eps contraction | Algebra | Levi-Civita contractions | Medium |
| PolarizationSum | Evaluate | Photon/gluon spin averaging | Easy |
| Phase space integration | Evaluate | Numerical cross-sections | Medium |
| Base.show methods | All | Interactive ergonomics | Easy |
| Symbolic N for colour | Algebra | Large-N, SU(2) comparison | Medium — NPoly type |

---

## Lessons Learned (READ THIS)

### The coefficient type IS the architecture
v1's `coeff::Any` + Expr building infected every function with adaptation code (~150 LOC of
`_mul_coeff`, `_flatten_product`, `_to_algsum`). v2's `DimPoly` eliminates all of it. Get the
coefficient representation right FIRST, before writing anything else. For new symbolic
quantities (N for colour, masses, etc.), define a proper algebraic type — never use `Any` or `Expr`.

### Separate fermion lines produce separate traces
`|M|² = Tr₁[...] × Tr₂[...]` is a PRODUCT of two AlgSums, not a single trace of a long
chain. The conjugate amplitude reverses gamma order and relabels indices (μ → μ_). This
was a subtle bug that gave wrong numerical answers.

### AlgSum for scalars, DiracExpr for matrices
"Everything returns AlgSum" was too aggressive. The self-energy numerator
γ^μ(m+γ·q)γ_μ = Dm·I + (2-D)·γ·q is matrix-valued — it can't live in AlgSum. DiracExpr
is `Vector{Tuple{AlgSum, DiracChain}}` — scalar coefficient × Dirac matrix. The two types
compose cleanly: `dirac_trace(::DiracExpr) → AlgSum`.

### Expand MomSumSlot before tracing
Don't add MomSumSlot dispatch to `gamma_pair`. Instead, expand `γ·(p-k)` into
`γ·p - γ·k` linearly before the trace algorithm runs. Keep the core algorithm simple.

### Research Julia idioms before each new layer
The Model layer initially used plain OOP structs. After researching Julia patterns, it
became: abstract interface (MultivariatePolynomials style), Field{Species} with traits,
gauge groups as types, callable FeynmanRules. The research step costs 5 minutes and saves
hours of refactoring.

---

## What To Do Next — Recommended Order

### Path A: Complete the 1-loop vertical (Stage B)

Add Layer 5 (Integrals) to the existing pipeline. Target: e+e- → μ+μ- at 1-loop.

1. Port PaVe{N} types from v1 to v2 (adapt to DimPoly coefficients)
2. Port tensor decomposition (tdec, rank 0-2) from v1
3. Port PaVe B-function reduction from v1
4. Implement Li₂ in pure Julia (~20 LOC)
5. Implement A₀, B₀ numerical evaluation (~100 LOC)
6. Build 1-loop vertex correction diagram
7. Test: Schwinger correction δσ/σ = (3α)/(4π)(π²/3 - 1/2)

See `src/v2/VERTICAL_PLAN.md` for detailed plan.

### Path B: Broaden the tree-level pipeline

Add more processes to validate the Model→Rules→Diagrams layers:

1. Compton scattering: e+γ → e+γ (needs PolarizationSum, 2 diagrams)
2. Bhabha scattering: e+e- → e+e- (s+t channels, same fermion species)
3. Møller scattering: e-e- → e-e- (identical particles)

Each process tests a different aspect of diagram generation and evaluation.

### Path C: Harden the algebra layer

Port v1's tested algorithms to v2's type system:

1. DiracTrick n≥3 (Mertig/Boehm/Denner general formula)
2. Eps contraction (Levi-Civita identity)
3. Port MUnit test translations from v1 (systematic coverage)

---

## Reference Codebases

All in `refs/` (gitignored):
- `refs/FeynCalc/` — Primary porting oracle. 186k LOC Mathematica. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/FormCalc/` — FormCalc reference (NOT directly ported).
- `refs/LoopTools/` — Loop integral numerics: dual FF+Denner Fortran implementation.

Architecture reports from Session 1 in `refs/reports/`.

---

## Quick Commands

```bash
# v1 tests (616 tests, should always pass)
julia --project=. -e 'using Pkg; Pkg.test()'

# v2 tests (123 tests, run individually — separate Julia processes)
julia --project=. test/v2/test_ee_mumu_x.jl
julia --project=. test/v2/test_coeff.jl
julia --project=. test/v2/test_colour.jl
julia --project=. test/v2/test_self_energy.jl
julia --project=. test/v2/test_vertical.jl

# Run all v2 tests sequentially
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Branch
git branch  # should show experimental/rebuild-v2

# Beads (may need restore)
bd init --force --prefix feynfeld && bd backup restore
bd stats
```
