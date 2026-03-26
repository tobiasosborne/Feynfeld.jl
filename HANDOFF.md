# HANDOFF — 2026-03-26 (Session 2)

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
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Run agents sequentially only.
    - **Clarification**: Rule 10 is about Julia precompilation. Read-only research/design agents CAN run in parallel.

**NEVER modify TensorGR.jl without explicit permission.** It is an active separate
project with its own workflow and handoff protocol.

---

## Current State

- **Phase 0 COMPLETE** (5/5 tasks, epic closed)
- **Phase 1a COMPLETE** (7/7 tasks, epic closed)
- **Phase 1b COMPLETE** (9/9 tasks, epic closed)
- **Phase 1c IN PROGRESS** (0/4 tasks — SU(N) colour, research + designs done)
- **Full test suite: 430 tests, ALL PASS**
- Beads: `bd ready` for available work. Note: beads DB may need `bd init --force --prefix feynfeld && bd backup restore` on fresh session.

---

## What Was Done This Session

### 1. Closed Phase 0 (tasks d4m.2–d4m.5)

Verified all Phase 0 implementation was complete and closed remaining issues.

### 2. Phase 1a: Pair Type System Design (feynfeld-43a.1)

Ran full Rule 4 workflow (3 agents):
- **Research agent**: Deep-dive into FeynCalc Pair/Momentum/Contract/ScalarProduct/Eps internals
- **Design A agent**: FeynCalc-faithful universal Pair approach
- **Design B agent**: Julia-idiomatic multiple-dispatch approach

**Decision: Hybrid** — Universal Pair type (from Design A, FeynCalc-proven) + SPContext explicit context (from Design B, thread-safe) + fold/reduce Contract (from Design B, no UpValues needed).

### 3. Implementation: Momentum, Pair, SPContext, Contract

New files created:
- `src/algebra/momentum.jl` (~45 LOC): Momentum(name, dim) type
- `src/algebra/pair.jl` (~130 LOC): Universal Pair with Orderless ordering + BMHV
- `src/algebra/sp_context.jl` (~55 LOC): Copy-on-write scalar product storage
- `src/algebra/contract.jl` (~155 LOC): Worklist reducer for Lorentz index contraction

Modified files:
- `src/algebra/types.jl`: LorentzIndex now uses DimSlot (was Union{Symbol,Int}). MetricTensor and ScalarProduct structs removed (replaced by Pair).
- `src/Feynfeld.jl`: Updated includes

Test files:
- `test/algebra/test_momentum.jl`: 14 tests
- `test/algebra/test_pair.jl`: 33 tests
- `test/algebra/test_sp_context.jl`: 10 tests
- `test/algebra/test_contract.jl`: 46 tests (full BMHV coverage)
- `test/algebra/test_types.jl`: Updated for LorentzIndex changes

### 4. Review findings addressed

Two reviewer agents ran. Key findings fixed:
- Added `Base.hash` for LorentzIndex (consistency with Momentum)
- Type annotation on `dim_contract` result to reduce boxing
- SPContext docstring corrected (copy-on-write, not immutable)
- Export cleanup (PairArg exported from momentum.jl where defined)
- Added edge case tests (self-contraction, cross-type hash, evanescent pairs)
- Added full BMHV test matrix (all 27 dim combinations for traces, partial, and vector contractions)

---

## Key Decisions / Lessons

### Carried from Session 1
- **FeynRules port MUST follow reference implementation**: Second-quantization operator algebra, NOT functional differentiation.
- **NEVER modify TensorGR.jl without explicit permission**
- **LorentzIndex is Feynfeld-specific, not a TIndex alias**: BMHV dimension tagging at index level
- **FCI/FCE unnecessary in Julia**: Constructors ARE the internal form

### New This Session
- **Pair not exported** (conflicts with Base.Pair): Users use SP/FV/MT convenience or `Feynfeld.Pair`
- **Pair is the universal building block**: Replaces MetricTensor, ScalarProduct structs
- **SPContext is explicit, not global**: Thread-safe, test-isolated. Passed as `ctx` kwarg to `contract()`
- **Contract uses fold/reduce**: No Mathematica UpValues trick. Explicit worklist with index inventory.
- **LorentzIndex default dim is Dim4()**: Matches FeynCalc convention. Old default was `:D`.
- **Read-only agents can run in parallel**: Rule 10 applies to Julia precompilation only.

---

## Known Limitations / Deferred

1. **`_mul_coeff` produces unsimplified symbolic Exprs**: `:(D * 4)` instead of `:(4D)`. Must fix before Dirac traces. Not blocking for Phase 1a.
2. **No expression tree for coefficients**: `contract(2 * FV(...), ...)` not expressible. Need `FMul`/`FAdd` expression algebra.
3. **MomentumSum not yet implemented**: Needed for ExpandScalarProduct (feynfeld-43a.4).
4. **Eps/EpsContract not yet implemented**: feynfeld-43a.6.

---

## TODO Next Session

1. **Phase 1c: SU(N) colour algebra** (`feynfeld-sem` epic, 4 tasks)
   - Research done, two design proposals done
   - Start with `feynfeld-sem.1`: SU(N) type system
   - Then SUNSimplify, SUNTrace, validation
2. Phase 1d (PaVe) is unblocked after 1c

## Quick Commands

```bash
bd ready                       # available work
bd show feynfeld-sem           # Phase 1c epic (0/4)
julia --project=. -e 'using Pkg; Pkg.test()'  # full test suite (430 tests)
```
