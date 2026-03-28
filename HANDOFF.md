# HANDOFF — 2026-03-28 (End of Session 5)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, anti-hallucination citation pattern, Julia idiom cheatsheet
2. Read `Feynfeld_PRD.md` — vision, spiral methodology, MUnit coverage target
3. Read `src/v2/DESIGN.md` — v2 type system, anti-patterns, cockroaches found
4. Run `bd init --force --prefix feynfeld && bd backup restore` if beads not working
5. Run `bd ready` to see available work
6. Run `for f in test/v2/test_*.jl; do julia --project=. "$f"; done` to verify 187 tests pass

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **ANTI-HALLUCINATION: CITE EVERYTHING.** Every formula must cite: local file
   path + equation number + verbatim copy. No exceptions. STOP if source missing.
3. **JULIA IDIOMATIC ALL THE WAY.** Read the cheatsheet in CLAUDE.md FIRST.
4. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.

**NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT STATE

### PRD v0.2 (rewritten this session)

Big vision: Julia-native, agent-facing, full-stack physics suite. Lagrangian →
predictions. Eventually ALL of physics (SUSY, quantum gravity, etc.). First
spirals: QED/QCD/EW validated against every FeynCalc MUnit test (15,222
assertions, ~10,000 translatable).

Development methodology: **spiral** (geometric mean of vertical × horizontal).
Each process drives both pipeline depth and MUnit coverage breadth.

### v2 codebase (branch `experimental/rebuild-v2`)

**187 tests pass across 7 test files. 2,500+ source LOC across 23 files.**

| Layer | Status | Key files |
|-------|--------|-----------|
| 4: Algebra | Tree-level COMPLETE | coeff, types, pair, expr, contract, expand_sp, dirac*, spin_sum, colour* |
| 1: Model | Hard-coded QED | model.jl (Field{Species}, GaugeGroup, traits) |
| 2: Rules | Hard-coded QED | rules.jl (FeynmanRules callable, dispatch on species) |
| 3: Diagrams | Hard-coded e+e-→μ+μ- | diagrams.jl (ExternalLeg, FeynmanDiagram) |
| 5: Integrals | PaVe types + A₀/B₀ | pave.jl, pave_eval.jl (QuadGK), schwinger.jl |
| 6: Evaluate | Tree-level σ | cross_section.jl (Problem/Solve, Mandelstam) |

### v1 (FROZEN — will be deleted)

616 tests in `src/algebra/`, `src/integrals/`. Mathematica-shaped anti-patterns
(tagged unions, `Any` coefficients, `Expr` building). Use ONLY as algorithmic
reference for porting. Do not extend. Do not import patterns from.

### Beads

79 total issues. 55 closed, 24 open. `bd --version` = 0.61.0, Dolt backend.
Server auto-starts. Open issues mapped to spiral plan.

---

## WHAT WAS DONE IN SESSION 5

### 1. Layer 5 (Integrals) — partial tracer bullet

Added PaVe{N} parametric type, numerical A₀/B₀ evaluation, Schwinger correction.

**New files:**
- `src/v2/pave.jl` (68 LOC) — PaVe{N} struct, named constructors A0/B0/B1/C0/D0
- `src/v2/pave_eval.jl` (90 LOC) — `evaluate` via dispatch, B0 via QuadGK Feynman parameter
- `src/v2/schwinger.jl` (53 LOC) — analytical Schwinger correction, vacuum polarization
- `test/v2/test_pave.jl` (158 LOC) — 51 tests: types, A₀, B₀ special cases, symmetry, general
- `test/v2/test_schwinger.jl` (64 LOC) — 13 tests: correction formula, vacuum pol, PaVe spot checks

**Dependencies added:** PolyLog.jl v2.6.2 (li2), QuadGK.jl v2.11.2 (adaptive quadrature)

**What this is NOT:** A full 1-loop vertical. The Schwinger correction function is
a hard-coded analytical constant, not computed from loop diagrams. C₀ (triangle),
tensor decomposition, renormalisation, and IR cancellation are all missing.

### 2. PRD v0.2 — complete rewrite

- Big vision: agent-facing full-stack physics, eventually all of physics
- Spiral methodology (process × MUnit coverage)
- Anti-hallucination citation rule (every formula: local path + eq number + verbatim)
- Julia idiom cheatsheet (replicated in PRD, JULIA_PATTERNS.md, CLAUDE.md)
- FeynCalc MUnit test coverage target (~10,000)
- v1 deletion policy

### 3. CLAUDE.md — updated

- New rules (anti-hallucination citation, Julia idiomatic)
- Beads + Dolt tracking instructions
- Ground truth acquisition (arXiv → TIB VPN → playwright-cli pattern)
- Julia idiom cheatsheet (third copy for LLM attention)

### 4. Beads triage

Closed 15 stale v1-era issues. 24 remain, mapped to spiral plan:
- Spiral 1 (Compton): PolarizationSum, topology generation
- Spirals 2-3 (Bhabha, QCD): DiracTrick full, SUN, amplitude squaring
- Spirals 4-6 (1-loop): C₀, D₀, PaVeReduce, TID, renormalisation
- Spiral 9+ (BSM/ULDM): Lagrangian DSL, SM model, ULDM portal

---

## KEY LEARNINGS (Session 5)

### 1. Don't hand-roll numerics

First attempt at B₀ used hand-derived analytical formulas for each special case.
The equal-mass formula had a sign bug (+β instead of -β) and missing +2. The
general formula was completely wrong. Reviewer agent caught it.

**Fix:** Replaced all special cases (except both-massless and one-massless)
with QuadGK adaptive quadrature on the Feynman parameter integral. Correct by
construction. ~14 digits accuracy. No hand-derived formulas to get wrong.

**Lesson:** Use Julia packages (QuadGK, PolyLog) instead of reimplementing.
This is now Rule 8 in CLAUDE.md.

### 2. PaVe is standalone, not in AlgFactor

Initial plan tried to make PaVe the 7th member of the AlgFactor union. Research
agent correctly identified this as OOP pollution: PaVe is scalar-valued, doesn't
carry Lorentz indices, doesn't contract or expand. It doesn't belong with tensor
factors. PaVe{N} is a standalone type with `evaluate` via dispatch.

**Where PaVe should live in AlgSum is an open question** for the learnings phase
after completing spirals 1-3.

### 3. The Schwinger correction requires C₀

A genuine 1-loop vertical (computing δσ/σ from diagrams) needs:
- C₀ (scalar triangle) for the vertex correction form factor
- Tensor decomposition to extract PaVe from loop amplitude
- Renormalisation (counterterms, Ward identity Z₁=Z₂)
- IR cancellation (virtual + real soft emission)

This is spirals 4-6, not a quick tracer bullet.

### 4. FeynCalc has 15,222 MUnit tests

Explored the full test suite. ~10,000 are translatable to Julia (the rest test
Mathematica-specific FCI/FCE conversions and typesetting). The tests encode
conventions, edge cases, and 25 years of cross-validation. They are the single
most valuable resource for the project.

### 5. Spiral methodology validated

The v2 rebuild proved vertical slices are more productive than horizontal layer
completion. The new spiral methodology combines both: each process is a vertical
spoke that drives horizontal MUnit coverage.

### 6. Anti-hallucination citation pattern is critical

The B₀ formula bug came from "I know the formula" without checking a reference.
Rule 2 now requires every formula to cite a local source with equation number
and verbatim copy. This is the most important rule for physics code.

---

## WHAT TO DO NEXT

### Spiral 1: Compton scattering (e+γ → e+γ)

This is the next process. It requires:

1. **PolarizationSum** — photon polarisation sum ε^μ ε^ν* → -g^{μν} (Feynman gauge)
   - MUnit tests: `refs/FeynCalc/Tests/Feynman/PolarizationSum.test`
   - Relatively simple, ~50 LOC

2. **Eps contraction** — Levi-Civita tensor contraction (if needed for cross-terms)
   - MUnit tests: `refs/FeynCalc/Tests/Lorentz/EpsContract.test`

3. **Two-diagram amplitude** — s-channel + u-channel
   - Extend diagrams.jl with Compton topology (or start proper topology generation)

4. **DiracTrick n≥3** — needed if the trace produces >2 sandwiched gammas
   - MUnit tests: `refs/FeynCalc/Tests/Dirac/DiracTrick.test` (577 tests!)

5. **Translate MUnit tests** for each function implemented

**Target:** `@test` validates Compton cross-section against Klein-Nishina formula.

### After Spiral 1

- Spiral 2: Bhabha (e+e-→e+e-) — t-channel, full DiracTrick
- Spiral 3: QCD (qq̄→gg) — full colour algebra
- Spiral 4: 1-loop vertex — C₀, tensor decomposition, PaVeReduce

### Ground truth papers still needed

The following papers should be acquired (TIB VPN) before spiral 4:
- Denner, Fortschr. Phys. 41 (1993) 307
- 't Hooft & Veltman, Nucl. Phys. B153 (1979) 365
- Passarino & Veltman, Nucl. Phys. B160 (1979) 151

P&S is already in `refs/papers/` (.djvu format).

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run all v2 tests (187 tests)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Run specific test
julia --project=. test/v2/test_pave.jl

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues

# Ground truth
ls refs/FeynCalc/Tests/         # MUnit test directories
ls refs/papers/                 # local paper copies
```
