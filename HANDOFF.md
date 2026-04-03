# HANDOFF — 2026-04-03 (Session 16: Spiral 10 Phases A-D, 1-loop box infrastructure)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 16 ACCOMPLISHMENTS

### 1. Spiral 10 Phases A-D: 1-loop box diagram infrastructure

Extended the 6-layer pipeline to support 1-loop box diagrams. Four new source
files + one test file, 551 LOC total. All files under 200 LOC limit.

**Phase A: LoopChannel types + box_channels()** (`loop_channels.jl`, 79 LOC)
- `LoopChannel` struct: topology, internal_fields, legs, loop_momentum
- `box_channels()`: enumerates direct + crossed box topologies for 2→2
- For ee→μμ QED: returns exactly 2 channels

**Phase B: Loop amplitude builder** (`loop_amplitude.jl`, 115 LOC)
- `build_loop_box_amplitude()`: constructs numerator DiracExpr chains
- `BoxDenominators`: stores accumulated momenta and PaVe invariant labels
- Momentum routing: D₀=q², D₁=(q+p₁)², D₂=(q+p₁+p₂)², D₃=(q+k₁_or_k₂)²
- Direct box: D₀(0,0,0,0,s,t,0,0,0,0), Crossed: D₀(0,0,0,0,s,u,0,0,0,0)

**Phase C: Tree × loop interference** (`loop_interference.jl`, 63 LOC)
- `spin_sum_tree_loop_interference()`: Σ_spins M_tree* × M_box
- Uses existing `_cross_line_trace()` from interference.jl
- Produces AlgSum with 33 terms containing SP(q,p_i), SP(q,q), external SPs
- After trace + contract + expand_sp: 6 rank-0, 21 rank-1, 6 rank-2 terms

**Phase D: Tensor Integral Decomposition (TID)** (`tid.jl`, 141 LOC)
- `evaluate_box_integral()`: evaluates ∫ d^D q N(q)/[D₀D₁D₂D₃]
- Rank 0: coefficient × D₀ (6 terms) ✓
- Rank 1: Σ_i (K_i·p_a) D_i using PV decomposition (21 terms) ✓
- Rank 2, SP(q,q): C₀⁽⁰⁾ (1 term) ✓
- Rank 2, SP(q,pa)×SP(q,pb): DEFERRED (5 terms, needs IR-regulated C-tensors)
- At √s=100 GeV, cosθ=0.5: finite ComplexF64 result (28/33 terms evaluated)

### 2. Known limitation: IR-divergent rank-2 TID

5 of 33 interference terms are rank-2 with SP(q,pa)×SP(q,pb) structure.
These require sub-triangle C₁/C₂ tensor evaluation which hits B₀(0,0,0)
(IR divergent in the massless case). Options to fix:
- (a) Use massive fermion regulation (m_e, m_μ > 0)
- (b) Use COLLIER for C-tensor evaluation (handles dim-reg)
- (c) Implement D₀₀/D_{ij} rank-2 tensor coefficients directly

Tracked in: feynfeld-544 (TID Phase D completion)

### 3. Beads

- **Claimed:** feynfeld-7h8 (1-loop amplitude builder) — Phases A-C DONE
- **Created:** feynfeld-544 (TID Phase D: rank-2 completion)
- **Created:** feynfeld-73g (Phase E+F: NLO evaluation + FeynCalc validation)
- **Dependency chain:** feynfeld-7h8 → feynfeld-544 → feynfeld-73g

---

## KNOWN ISSUES AND BLOCKERS

### P1: Rank-2 TID IR divergence (feynfeld-544)
5 rank-2 terms skipped in massless box. Blocks full NLO validation.

### P1 (pre-existing): Δα imaginary part — 1 flaky test
Intermittent. Not related to Session 16 changes.

### P3: `propagator_num(::Boson)` is dead code
Pre-existing. Defined but never called.

### P3: `PropagatorRule` exported but undefined
Pre-existing.

---

## WHAT TO DO NEXT

### Priority 1: Complete rank-2 TID (feynfeld-544)
Fix the IR-divergent sub-triangle C₁/C₂ in the massless box.
Best approach: use COLLIER for C-tensor evaluation (handles dim-reg natively).
Alternative: add small fermion masses for IR regulation.

### Priority 2: NLO validation (feynfeld-73g)
Once TID is complete, evaluate the full box interference at multiple (s,t)
points and compare against FeynCalc ElAel-MuAmu.m reference.
Ref: arXiv:hep-ph/0010075, Eq. 2.32.

### Priority 3: MUnit test porting (~370 remaining)
Continues alongside spirals. High-value targets:
- feynfeld-rcd: Contract D-dim section (20 tests)
- feynfeld-4mm: SUNSimplify (78 tests)
- feynfeld-36h: SUNTrace (24 tests)

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
8. **TIERED WORKFLOW.** Core (>20 LOC): 3 research + 1 review. Small: 1+1. Trivial: direct.
9. **NEVER modify core algebra files without explicit permission.**

---

## PROJECT STATE

### Branch and code location
- **Branch:** `master`
- **v2 source:** `src/v2/` (43 files, ~4,600 LOC)
- **v2 tests:** `test/v2/` (21 files + munit/) — 540+ tests
- **v1:** FROZEN. Do not extend.

### New/modified files in Session 16

| File | LOC | What |
|------|-----|------|
| `src/v2/loop_channels.jl` | 79 | NEW: LoopChannel type + box_channels() |
| `src/v2/loop_amplitude.jl` | 115 | NEW: build_loop_box_amplitude() + BoxDenominators |
| `src/v2/loop_interference.jl` | 63 | NEW: spin_sum_tree_loop_interference() |
| `src/v2/tid.jl` | 141 | NEW: evaluate_box_integral() + TID |
| `src/v2/FeynfeldX.jl` | ~130 | MODIFIED: includes + exports for loop infrastructure |
| `test/v2/test_box_ee_mumu.jl` | 153 | NEW: 136 tests for Phases A-D |

### Box diagram momentum routing

For the direct box e⁻(p₁)+e⁺(p₂)→μ⁻(k₁)+μ⁺(k₂):
- D₀ = q² (photon₁)
- D₁ = (q+p₁)² (electron internal)
- D₂ = (q+p₁+p₂)² (photon₂)
- D₃ = (q+k₁)² (muon internal)
- PaVe: D₀(0,0,0,0,s,t, 0,0,0,0)

For the crossed box: D₃ = (q+k₂)², PaVe: D₀(0,0,0,0,s,u, 0,0,0,0)

### Accumulated momenta (PaVe convention)
- K₁ = p₁
- K₂ = p₁ + p₂  (MomentumSum)
- K₃ = k₁ (direct) or k₂ (crossed)

---

## QUICK COMMANDS

```bash
# Run single test (fast)
julia --project=. test/v2/test_box_ee_mumu.jl   # 5s, Spiral 10 infrastructure

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Beads
bd ready              # available work
bd show feynfeld-544  # TID completion
bd show feynfeld-73g  # NLO validation
bd stats              # project health

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
