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
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Read-only agents CAN run parallel.

**NEVER modify TensorGR.jl without explicit permission.**

---

## Current State

- **Phase 0 COMPLETE** (5/5 tasks, epic closed)
- **Phase 1a COMPLETE** (7/7 tasks, epic closed) — Lorentz: Pair, Contract, ExpandSP, Eps
- **Phase 1b COMPLETE** (9/9 tasks, epic closed) — Dirac: DiracGamma, DiracTrick, DiracTrace, DiracSimplify
- **Phase 1c COMPLETE** (4/4 tasks, epic closed) — SU(N): SUNT, SUNF, SUND, SUNTrace, SUNSimplify
- **Phase 1d COMPLETE** (6/6 tasks, epic closed) — PaVe: PaVe{N}, FAD, Tdec, PaVeReduce
- **Full test suite: 506 tests, ALL PASS**
- **Closed: 45/77 beads issues**

---

## What Was Done This Session

### Session 2 completed ALL of Phase 1 (algebra layer):

1. **Phase 1a (Lorentz)**: Pair universal type, Momentum/MomentumSum, Contract (full BMHV 27-dim matrix), ExpandScalarProduct, Eps/EpsContract, SPContext. 7/7 tasks.

2. **Phase 1b (Dirac)**: DiracGamma tagged slots (LI/Mom/Special), DiracChain/DOT, DiracTrick (g5², projectors, slash², trace, sandwich n=1,2), DiracTrace recursive formula + chiral g5 base case, DiracOrder, DiracEquation, DiracSimplify orchestrator, DiracScheme (NDR/BMHV/Larin enum). 9/9 tasks.

3. **Phase 1c (SU(N))**: SUNIndex/SUNFIndex, SUNT, SUNTF, SUNF (antisymmetric auto-sorted), SUND (symmetric), SUNDelta/SUNFDelta, ColourChain, SUNTrace recursive, delta_trace, contract_ff/dd/fd. 4/4 tasks.

4. **Phase 1d (PaVe)**: PaVe{N} parametric, A0/A00/B0/B1/B00/B11/C0/D0, FeynAmpDenominator + 3 propagator types, Tdec rank 0/1/2, PaVeReduce B-functions. 6/6 tasks.

### Source files created this session (all under 200 LOC):

| File | LOC | Purpose |
|------|-----|---------|
| `src/algebra/momentum.jl` | 130 | Momentum, MomentumSum, arithmetic |
| `src/algebra/pair.jl` | 135 | Universal Pair, BMHV, SP/FV/MT |
| `src/algebra/sp_context.jl` | 55 | Copy-on-write SP storage |
| `src/algebra/contract.jl` | 155 | Lorentz index contraction |
| `src/algebra/expand_sp.jl` | 60 | Bilinear expansion |
| `src/algebra/eps.jl` | 135 | Levi-Civita + EpsContract |
| `src/algebra/dirac_types.jl` | 150 | DiracGamma, Spinor, DiracChain |
| `src/algebra/dirac_chain.jl` | 80 | DOT constructor, dot_simplify |
| `src/algebra/dirac_trick.jl` | 185 | Core Dirac simplification rules |
| `src/algebra/dirac_trace.jl` | 110 | Recursive trace formula |
| `src/algebra/dirac_order.jl` | 55 | Normal ordering |
| `src/algebra/dirac_equation.jl` | 75 | Dirac equation at boundaries |
| `src/algebra/dirac_simplify.jl` | 40 | Master orchestrator |
| `src/algebra/dirac_scheme.jl` | 45 | NDR/BMHV/Larin scheme |
| `src/algebra/colour_types.jl` | 130 | SU(N) types |
| `src/algebra/colour_trace.jl` | 75 | SUNTrace recursive |
| `src/algebra/colour_simplify.jl` | 100 | Structure constant contractions |
| `src/integrals/pave.jl` | 85 | PaVe{N} + named constructors |
| `src/integrals/feynamp_denominator.jl` | 90 | FAD + propagator types |
| `src/integrals/tdec.jl` | 75 | Tensor decomposition rank 0/1/2 |
| `src/integrals/pave_reduce.jl` | 55 | B-function PaVe reduction |

---

## Known Limitations / Deferred

1. **Symbolic coefficient simplification**: `_mul_coeff` produces unsimplified `Expr` trees
2. **No expression tree for coefficients**: `contract(n * FV(...))` not expressible
3. **DiracTrick sandwich n≥3**: General contraction formula not implemented
4. **DiracTrace chiral n>4**: Recursive chiral trace for 6+ gammas with g5
5. **BMHV/Larin scheme**: Registered but not fully wired into DiracTrick
6. **SUNSimplify**: Cvitanovic/Fierz completeness relation not yet implemented
7. **Tdec rank ≥3**: Higher-rank tensor decomposition
8. **PaVeReduce C/D**: Only B-function reductions implemented
9. **ToPaVe**: FAD → PaVe conversion not yet implemented

---

## TODO Next Session

1. **Phase 2**: LoopTools numerical integration (`feynfeld-62k` epic)
2. **Phase 3**: Model + Rules — FeynRules port (`feynfeld-ntj` epic)
3. **Phase 4**: Diagrams — FeynArts port (`feynfeld-c0n` epic)
4. **Phase 5**: Evaluate — amplitude squaring (`feynfeld-89y` epic)
5. **Phase 6**: ULDM application (`feynfeld-42d` epic)

## Quick Commands

```bash
bd stats                       # project statistics
bd ready                       # available work
julia --project=. -e 'using Pkg; Pkg.test()'  # 506 tests
```
