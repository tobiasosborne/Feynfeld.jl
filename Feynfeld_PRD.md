# Feynfeld.jl — Product Requirements Document

**Version:** 0.2
**Date:** 2026-03-28
**Author:** Tobias J. Osborne / LUH
**Status:** Active development (v2 on `experimental/rebuild-v2`)

---

## 1. Vision

Feynfeld.jl is a Julia-native, agent-facing, full-stack physics computation suite.

You start Claude Code with a physics idea. You discuss the Lagrangian. Then everything
is automated: the result is predictions of physical observables. Full stack, from
fundamental physics idea to prediction. No Mathematica, no Fortran FFI, no manual
stitching of notebooks.

```
Human: "I want the one-loop correction to e+e- → μ+μ- in QED."

Agent: reads Lagrangian → extracts Feynman rules → generates diagrams →
       contracts indices → traces gamma matrices → reduces PaVe integrals →
       evaluates numerically → returns σ_NLO with uncertainty.
```

### 1.1 Scope trajectory (updated 2026-03-29)

**Done (spirals 0–7):** QED tree-level (e+e-→μ+μ-, Compton, Bhabha), QCD tree
(qq̄→gg), EW tree (e+e-→W+W-), one-loop vertex/self-energy/vacuum polarization.
301 tests, 28 source files.

**Current (spiral 8):** Bug fixes, γ5 traces, Eps contraction. Unblock chiral
physics and the full EW sector.

**Next (spirals 9–10):** Pipeline consolidation. Every existing process through
the full 6-layer pipeline. Minimal diagram generation for tree-level 2→2 via
channel filtering (~200 LOC). EW/QCD models with explicit SM vertex dispatch.
D₀ evaluation and box diagrams.

**Medium-term (spirals 11–15):** Full SM at one loop. Automated 1-loop tensor
reduction (TID/OPP). Renormalization. 2→3 tree-level diagrams. IR subtraction
(Catani-Seymour or FKS). PDF convolution interface. Phase space Monte Carlo.

**Long-term (spirals 16+):** BSM via Lagrangian DSL / UFO import. ULDM portal
for VLBAI. 2-loop via IBP reduction. TensorGR.jl bridge for graviton exchange.

### 1.2 The endgame interaction

```
Human: "Compute pp → H + jet at NLO QCD."

Agent: parses process → identifies partonic subprocesses (gg→Hg, qg→Hq, qq̄→Hg)
     → generates tree + 1-loop diagrams for each
     → contracts, traces, reduces to PaVe
     → evaluates loop integrals numerically
     → adds real emission with IR subtraction
     → convolves with PDFs
     → integrates over phase space
     → returns σ(pp → H+jet) = 29.6 ± 0.3 pb at √s = 13 TeV
```

This is the target. It requires every layer of the pipeline to work for
arbitrary SM processes at NLO. The architectural foundation (pipeline
principle, Julia dispatch, DimPoly coefficients, Dict-based AlgSum) scales
to this level. The limiting factors are physics algorithms (IBP reduction,
IR subtraction), not software architecture.

**Hardness scale for loop orders:**

| Order | Status | Effort | Notes |
|-------|--------|--------|-------|
| Tree (LO) | Spiral 9 | Weeks | Channel filtering, ~200 LOC |
| 1-loop (NLO) | Spirals 10-14 | Months | PV reduction + IR subtraction |
| 2-loop (NNLO) | Spirals 15+ | Person-years | IBP reduction (FIRE/Kira-class) |
| 3-loop (N3LO) | Aspirational | Years | Cutting-edge research, handful of results exist |
| 4-loop+ | Out of scope | — | Nobody has done this for physical collider processes |

### 1.3 Agent-facing, human-convenient

The primary user is an AI agent (Claude Code) that reads the Lagrangian, plans
the calculation, and calls Feynfeld.jl functions. The API is designed for
programmatic composition, not interactive exploration.

But humans work at the REPL too. Convenience constructors (`SP(:p,:q)`,
`GA(:mu)`, `qed_model()`) and `Base.show` methods make interactive use natural.

### 1.3 Relationship to TensorGR.jl

TensorGR.jl (Julia port of xAct for tensor calculus in GR) is a sibling
project. They share no code at present — v2 has its own Lorentz algebra (~200
LOC). Future bridge: graviton exchange, scalar fields on curved backgrounds,
quantum gravity calculations. The bridge is a future spiral, not a dependency.

---

## 2. Architecture

### 2.1 One package, six layers

```
Layer 1: Model      Lagrangian, fields, parameters, gauge groups
Layer 2: Rules      Lagrangian → Feynman rules (second quantisation)
Layer 3: Diagrams   Topology generation, field insertion, amplitude construction
Layer 4: Algebra    Lorentz · Dirac · Colour · Tensor contraction/traces
Layer 5: Integrals  PaVe scalar functions, tensor decomposition, numerical eval
Layer 6: Evaluate   |M|², cross-sections, decay rates, RGE, observables
```

**The Algebra layer is the type system.** Everything above constructs expressions
in its types. Everything below consumes them.

**THE PIPELINE PRINCIPLE: The pipeline IS the architecture.** If a process
bypasses it, the architecture is incomplete. Every physics process must flow
through all 6 layers. Standalone analytical recipes are reference
implementations for cross-validation, not substitutes. No hand-built
amplitudes in test files. Fix the pipeline, don't work around it.

### 2.2 The coefficient type IS the architecture

The single most important design insight from v2: `DimPoly` (polynomial in D)
as the coefficient type eliminated ~150 LOC of adaptation code from v1. Get the
coefficient representation right FIRST.

```
Coeff = Union{Rational{Int}, DimPoly}    # never Any, never Expr
AlgSum = Dict{FactorKey, Coeff}          # O(1) like-term collection
```

For new symbolic quantities (N for colour, masses, couplings), define a proper
algebraic type. Never use `Any` or `Expr`.

### 2.3 What exists (v2, 301 tests, updated 2026-03-29)

| Layer | Status | Files |
|-------|--------|-------|
| 4: Algebra | Tree+1-loop complete | coeff, types, pair, expr, contract, expand_sp, dirac*, spin_sum, colour*, polarization_sum |
| 1: Model | Hard-coded QED | model.jl |
| 2: Rules | Hard-coded QED | rules.jl |
| 3: Diagrams | Hard-coded e+e-→μ+μ- (scaffolding) | diagrams.jl |
| 5: Integrals | PaVe A₀/B₀/B₁/C₀/C₁/C₂ + standalone recipes | pave.jl, pave_eval.jl, schwinger.jl, vertex.jl, running_alpha.jl |
| 6: Evaluate | Tree-level σ + EW cross-section | cross_section.jl, ew_parameters.jl, ew_cross_section.jl |

28 source files (~3,400 LOC), 16 test files (~2,000 LOC), 301 tests.
All code in `src/v2/`. Design choices documented in `src/v2/DESIGN.md`.

**Architectural review note:** Layer 5 standalone recipes (schwinger.jl, vertex.jl,
running_alpha.jl, ew_cross_section.jl) are analytical reference implementations,
not pipeline components. They validate physics but do not exercise the six-layer
pipeline. Future spirals should add pipeline tests that reproduce these results
through Model→Evaluate.

### 2.4 v1 is archived

v1 (616 tests, `src/algebra/`, `src/integrals/`) is **frozen and will be
deleted** once v2 achieves equivalent coverage. v1 mirrored FeynCalc's
Mathematica patterns (tagged unions, `Expr` coefficients, polymorphic returns)
— these are anti-patterns in Julia. v1 exists only as an algorithmic reference
for porting. Do not extend. Do not import patterns from.

---

## 3. Development Methodology: The Spiral

### 3.1 Geometric mean of vertical and horizontal

Neither pure vertical (one process end-to-end) nor pure horizontal (one layer
complete) is optimal. The spiral combines both:

1. **Pick a process** (vertical target, e.g., Compton scattering)
2. **Identify the FeynCalc functions** it requires
3. **Translate the MUnit tests** for those functions (horizontal coverage)
4. **Implement in v2** until the MUnit tests pass
5. **Run the process end-to-end** (vertical validation)
6. **Measure coverage gap** → pick next process to maximise coverage

Each process is a spoke that drives outward through the MUnit test suite. After
enough spokes, you've covered the full algebra.

### 3.2 FeynCalc MUnit tests as ground truth

FeynCalc has 15,222 MUnit test assertions across 300 files. ~10,000 are
directly translatable to Julia (the rest test Mathematica-specific infrastructure).

The MUnit suite solves the three hardest problems:
- **Convention choices**: metric signature, gamma trace normalisation, PaVe
  argument ordering — all encoded implicitly in the tests.
- **Edge cases**: zero masses, collinear momenta, degenerate limits.
- **25 years of cross-validation**: the FeynCalc team has validated against
  textbooks, LoopTools, FORM, and published papers.

Every translated MUnit test includes a comment citing the source file and test ID:
```julia
# Source: refs/FeynCalc/Tests/Lorentz/Contract.test, Test #42
@test contract(MT(:mu,:nu) * MT(:mu,:rho)) == MT(:nu,:rho)  # fcstContract-ID42
```

### 3.3 Spiral sequence (updated 2026-03-29 after architectural review)

| Spiral | Process | Status | Tests |
|--------|---------|--------|-------|
| 0 | e+e-→μ+μ- tree | DONE (Session 4) | ~65 |
| 1 | Compton e+γ→e+γ | DONE (Session 5) | +58 |
| 2 | Bhabha e+e-→e+e- | DONE (Session 6) | +4 |
| 3 | QCD qq̄→gg | DONE (Session 6) | +2 |
| 4 | 1-loop self-energy Σ(p) | DONE (Session 6) | +13 |
| 5 | 1-loop vertex (g-2) | DONE (Session 7) | +32 |
| 6 | Running α(q²) | DONE (Session 7) | +36 |
| 7 | EW e+e-→W+W- tree | DONE (Session 7) | +33 |
| **8** | **Bug fixes + γ5 + Eps** | **NEXT** | target: +50 |
| 9 | Diagram generation + EW model | Planned | target: +50 |
| 10 | D₀ + box diagrams | Planned | target: +50 |
| 11+ | BSM / ULDM | Planned | — |

**Revised from Session 8 review:** Spiral 8 was "MUnit mop-up" — changed to
"bug fixes + γ5 + Eps contraction" based on 6-agent architectural review.
MUnit translation is now continuous alongside spirals, not a dedicated phase.

**Metric change:** Drop raw test count target (was ~10,000). Adopt function
coverage: track how many of ~60 core FeynCalc functions have ≥5 translated
MUnit tests. Current: ~6 functions with ≥5 tests (~10% function coverage).

---

## 4. Rules

These are non-negotiable. They are learnt through suffering.

### TOBIAS'S RULES — FOLLOW TO THE LETTER

1. **GROUND TRUTH = PHYSICS.** Not pinned numbers. Not LLM memory. Not "I think
   the formula is." Physics is the ONLY truth.

2. **ANTI-HALLUCINATION: CITE EVERYTHING.** Every physics formula in source code
   AND test code must have a comment citing: (a) the local file path of the
   reference, (b) the exact equation number, (c) a verbatim copy of the equation
   as it appears in the source. No equation is valid without this triple.
   ```julia
   # Ref: refs/papers/Denner1993.pdf, Eq. (4.18)
   # "B₁(p², m₀², m₁²) = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)"
   B1 = (a0_m1 - a0_m0 - (p2 + m02 - m12) * b0) / (2 * p2)
   ```
   If you cannot find the formula in a local source, STOP and ask for it. Do not
   proceed with an uncited formula. This is the critical anti-hallucination pattern.

3. **SKEPTICISM.** All subagent work, handoffs — verify everything twice.

4. **ALL BUGS ARE DEEP.** Deep, complex, interlocked. Do not underestimate. Do
   not apply bandaids. Best-practices full solutions only.

5. **WORKFLOW.** 3 subagents before any core code change: (a) research the local
   copy of the porting source, (b) solution proposal 1, (c) solution proposal 2.
   Then choose the better approach.

6. **REVIEW.** Rigorous reviewer agent after every core change. No exceptions.

7. **TESTING.** Targeted only, or full suite in background. Every exported
   function has at least one `@test`. MUnit test translation protocol (§3.2).

8. **JULIA IDIOMATIC ALL THE WAY.** Read the Julia Idiom Cheatsheet (§6) before
   writing any code. If you catch yourself writing `x.field isa SomeType`, stop.
   Use existing Julia packages (QuadGK, PolyLog, etc.) — never reimplement
   standard library functionality.

9. **NO PARALLEL JULIA AGENTS.** Precompilation cache conflicts. Read-only
   research/design agents CAN run in parallel.

10. **RESEARCH IDIOMS FIRST.** Before every new layer, research Julia-idiomatic
    patterns for 5 minutes. Check how ecosystem packages solve the same problem.
    This saves hours of refactoring.

11. **LOC LIMIT.** No source file exceeds ~200 lines. Split aggressively.

12. **REPEAT RULES.** Repeat occasionally to maintain focus.

### Anti-hallucination reference infrastructure

All reference sources must be stored locally:
- `refs/FeynCalc/` — Primary porting oracle (MUnit tests in `Tests/`)
- `refs/FeynArts/` — Diagram generation reference
- `refs/FeynRules/` — Model/Lagrangian reference
- `refs/LoopTools/` — Loop integral numerics (Fortran source)
- `refs/papers/` — Published papers (Denner 1993, 't Hooft-Veltman 1979, P&S, etc.)
- `refs/reports/` — Architecture analysis reports

If a reference is not locally available, download it before citing.

---

## 5. Design Principles

### 5.1 Types over patterns

Mathematica uses rule-based pattern matching. Julia uses multiple dispatch on
algebraic types. The translation:

- `f[x_?SomeTest] := ...` → `f(x::SomeType) = ...`
- `expr /. {rule1, rule2}` → method dispatch or explicit `simplify` passes
- Tagged unions with `isa` checks → parametric types with dispatch

### 5.2 Immutable expressions, functional transformations

All algebraic expressions are immutable. Simplification produces new expressions,
never mutates. Canonical forms are enforced in constructors (sort, normalise).

### 5.3 Uniform return types

- Scalar algebra → `AlgSum` (always)
- Matrix-valued Dirac → `DiracExpr` (always)
- No mixed returns (`Int|Pair|Expr|Vector`)

### 5.4 Construction is cheap, simplification is explicit

Constructors canonicalise (sort, normalise signs) but never simplify. Contraction,
trace, expand are explicit function calls.

### 5.5 Use Julia packages

- **QuadGK.jl** for numerical integration (not hand-rolled quadrature)
- **PolyLog.jl** for dilogarithms (not reimplemented Li₂)
- **LoopTools.jl** for cross-validation (test-time only)
- **ScopedValues** for implicit context (SPContext pattern)

---

## 6. Julia Idiom Cheatsheet

**READ THIS BEFORE WRITING ANY CODE.** This section is replicated in
`JULIA_PATTERNS.md` and `CLAUDE.md` to ensure every agent sees it.

### DO

```julia
# Parametric types for dispatch (not tagged unions)
struct DiracGamma{S<:DiracSlot}
    slot::S
end

# Multiple dispatch (not if/elseif isa cascades)
contract(a::MetricTensor, b::FourVector, idx) = ...
contract(a::FourVector, b::FourVector, idx) = ...

# Dict-based expression storage (O(1) like-term collection)
struct AlgSum
    terms::Dict{FactorKey, Coeff}
end

# Concrete small Unions (Julia optimises up to ~4 types)
const Coeff = Union{Rational{Int}, DimPoly}

# Named constructors (plain functions returning types)
A0(m2::Real) = PaVe{1}(Int[], Float64[], [Float64(m2)])

# ScopedValues for implicit context
const CURRENT_SP = ScopedValue(SPContext())

# Use existing packages
using QuadGK: quadgk
using PolyLog: li2

# evaluate via dispatch, not type-in-function-name
evaluate(pv::PaVe{1}; mu2=1.0) = ...
evaluate(pv::PaVe{2}; mu2=1.0) = ...
```

### DO NOT

```julia
# ✗ Tagged union with isa checks
if x isa ScalarProduct ... elseif x isa MetricTensor ... end

# ✗ Any or Expr as coefficient type
struct BadTerm; coeff::Any; end

# ✗ Building Expr trees
coeff = :($(a) * $(b) + $(c))

# ✗ Type-in-function-name (Java pattern)
evaluate_pave(x)     # ✗
evaluate(x::PaVe)    # ✓

# ✗ Wrapper structs when a function will do
struct SchwingerProblem; alpha::Float64; end   # ✗
schwinger_correction(; alpha=1/137.036) = ...  # ✓

# ✗ Hand-rolled numerical methods
const MY_GAUSS_NODES = (...)  # ✗
quadgk(f, 0, 1)              # ✓

# ✗ Reimplementing standard functions
function my_dilog(z) ... end  # ✗
using PolyLog: li2            # ✓

# ✗ OOP struct hierarchies by default
abstract type AbstractIntegral end  # ✗ (unless dispatch genuinely needed)

# ✗ Forcing types into unions they don't belong in
const AlgFactor = Union{Pair, Eps, ..., PaVe}  # ✗ PaVe is scalar, not a tensor factor

# ✗ Global mutable state
CURRENT_MODEL = nothing  # ✗
const CURRENT_MODEL = ScopedValue(nothing)  # ✓
```

### Julia ecosystem patterns to follow

| Pattern | Example package | What to learn |
|---------|----------------|---------------|
| Dict-based expression storage | SymbolicUtils.jl | `AlgSum` pattern |
| Problem/Solve | DiffEq.jl | `CrossSectionProblem` → `solve_tree` |
| Abstract interface + traits | MultivariatePolynomials.jl | `AbstractModel` pattern |
| Callable structs | Flux.jl | `FeynmanRules(fields)` |
| ScopedValues | Julia 1.11+ stdlib | `SPContext` pattern |

---

## 7. Success Criteria

### 7.1 Minimum viable (DONE ✓)

e+e- → μ+μ- tree-level from Model → Rules → Diagrams → Algebra → Evaluate.
P&S Eq. (5.10) and (5.12) to machine precision. **Achieved in v2.**

### 7.2 Spiral 4 milestone

1-loop vertex correction computed FROM the pipeline (not hard-coded). PaVe
integrals extracted, reduced, evaluated. ~5,000 MUnit tests passing.

### 7.3 ULDM milestone

Scalar ULDM portal Lagrangian → coupling constants d_{m_e}, d_e → one-loop
RGE → numerical values for VLBAI. Connect to downstream Julia pipeline.

### 7.4 Full SM milestone

All 2→2 SM processes at tree level and one loop. ~10,000 MUnit tests passing.
All FeynCalc physics tests covered.

### 7.5 Community milestone

Reproduce AION sensitivity projections from Badurina et al. (2109.10965).

---

## 8. Testing Strategy

### 8.1 Ground truth hierarchy

1. **Published papers and textbooks** (locally stored, cited by equation number)
2. **FeynCalc MUnit tests** (primary porting oracle, ~10,000 translatable)
3. **Cross-validation** against LoopTools, COLLIER, or Package-X (≥2 checks per
   numerical function)
4. **Physics invariants**: Ward identities, gauge invariance, unitarity cuts

Rule 6 applies absolutely: if MUnit and textbook disagree, textbook wins.

### 8.2 MUnit translation protocol

For each FeynCalc MUnit test:
1. Read the `.test` file from `refs/FeynCalc/Tests/`
2. Translate each `Test[]` to a Julia `@test`, preserving math exactly
3. Document the source file and test ID in a comment
4. If MUnit test and textbook disagree, the textbook wins (Rule 1)
5. Cite the textbook equation that validates the test

### 8.3 Anti-hallucination test pattern

Every `@test` that checks a physics result must cite its source:
```julia
# Source: refs/FeynCalc/Tests/Lorentz/Contract.test, Test #42
# Cross-ref: P&S Eq. (A.27): g^μν g_μρ = δ^ν_ρ
@test contract(MT(:mu,:nu) * MT(:mu,:rho)) == ...
```

---

## 9. Risks

| Risk | Severity | Mitigation |
|------|----------|------------|
| Dirac algebra canonicalisation subtler than expected | High | MUnit tests catch it. Anti-panic: 3 subagents before fix. |
| FeynCalc MUnit tests contain bugs | Medium | Rule 1. Cross-validate against textbooks. |
| Convention drift between layers | High | MUnit tests enforce conventions end-to-end. |
| Agent hallucinating physics formulas | Critical | Rule 2. Every formula cited from local source. |
| Julia package ecosystem changes | Low | Pin versions in Project.toml. |
| 1-loop vertical genuinely hard | High | Spirals 4–6 are dedicated to it. Don't rush. |
| v1 patterns leaking into v2 | Medium | v1 archived. DESIGN.md anti-patterns list. |

---

## 10. References

### Primary sources (must be locally available)

- M. E. Peskin, D. V. Schroeder, *An Introduction to QFT* (1995)
- A. Denner, Fortschr. Phys. 41 (1993) 307 — one-loop techniques
- G. 't Hooft, M. Veltman, Nucl. Phys. B153 (1979) 365 — scalar integrals
- G. Passarino, M. Veltman, Nucl. Phys. B160 (1979) 151 — PV reduction
- V. Shtabovenko et al., "FeynCalc 10", arXiv:2312.14089

### Reference codebases (in `refs/`, gitignored)

- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation
- `refs/FeynRules/` — Model/Lagrangian
- `refs/LoopTools/` — Loop integral numerics (Fortran)
- `refs/papers/` — Local copies of cited papers

### ULDM-specific

- Badurina et al., Phys. Rev. D 105 (2022) 023006, arXiv:2109.10965
- Badurina et al., arXiv:2108.02468
- Badurina et al., arXiv:1911.11755

---

## Appendix A: Name

*Feynfeld* = Feynman + Feld (German: field).
