# HANDOFF — 2026-04-03 (Session 16: Spiral 10 complete, Δα fix, TID rank-2 via COLLIER)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (406 tests, ~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 16 ACCOMPLISHMENTS

### 1. Spiral 10 Phases A-D: 1-loop box diagram infrastructure (COMPLETE)

Extended the 6-layer pipeline to support 1-loop box diagrams. Five new source
files + one test file. All files under 200 LOC limit.

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
- Produces AlgSum with 33 terms: 6 rank-0, 21 rank-1, 6 rank-2

**Phase D: Tensor Integral Decomposition** (`tid.jl`, ~170 LOC)
- `evaluate_box_integral()`: evaluates ∫ d^D q N(q)/[D₀D₁D₂D₃]
- Rank 0: coefficient × D₀ via COLLIER ✓
- Rank 1: Σ_i (K_i·p_a) D_i via PV decomposition ✓
- Rank 2, SP(q,q): C₀⁽⁰⁾ + m₀²D₀ ✓
- Rank 2, SP(q,pa)×SP(q,pb): PV cancellation → sub-triangle C₁/C₂ ✓
- **ALL 33/33 terms evaluate.** Finite ComplexF64 at all kinematics.

### 2. Rank-2 TID IR fix: COLLIER B₀ fallback

The rank-2 TID was blocked by B₀(0,0,0) — IR-divergent in the massless box.
Sub-triangle C₁/C₂ PV reduction requires B₀ calls that can hit (0,0,0).

**Fix** (5 LOC in `pave_eval.jl`): when `B₀(0,0,0)` is requested, fall back to
COLLIER's `b0_coli_` which handles dimensional regularization natively and
returns the finite MS-bar remainder (= 0 for this degenerate case).

The `b0_coli_` symbol is a plain (non-module) Fortran function in libcollier.so,
callable directly via ccall. No Fortran module name mangling needed.

### 3. Δα imaginary part sign fix

The test `@test imag(da_mz) > 0.0` was **wrong**, not flaky.
Physics: Δα = -Π̂ (Denner Eq. 3.10). Absorptive Im(Π̂) > 0 → Im(Δα) < 0.
Fixed to `@test imag(da_mz) < 0.0`. All 406 tests now pass, zero flaky.

### 4. Beads closed

- **feynfeld-7h8** CLOSED: 1-loop amplitude builder (Phases A-C)
- **feynfeld-544** CLOSED: Rank-2 TID complete (COLLIER B₀ fallback)
- **feynfeld-8c4** CREATED: Epic for world-class diagram generation (P4 backlog)

---

## KNOWN ISSUES

### P3: `propagator_num(::Boson)` is dead code
Pre-existing. Defined but never called.

### P3: `PropagatorRule` exported but undefined
Pre-existing.

---

## WHAT TO DO NEXT

### Priority 1: NLO box validation (feynfeld-73g)

The full box integral pipeline works (all 33 terms, finite ComplexF64).
**Next step**: evaluate at multiple (s,t) points and cross-validate against
FeynCalc `ElAel-MuAmu.m` reference (arXiv:hep-ph/0010075, Eq. 2.32).

Needs:
1. Assemble coupling factors (e⁶ for tree×box interference)
2. Sum direct + crossed box contributions
3. Add tree propagator denominator (1/s from tree photon)
4. Evaluate 2·Re(M_tree* × M_box) at √s = 50, 100, 200 GeV
5. Compare box coefficient of D₀(s,t) against FeynCalc GLI entries
6. The FeynCalc reference separates into 3 topologies:
   - GLI["fctopology1", {1,1,1,1}] → box with (s,t)
   - GLI["fctopology2", {1,1,1,1}] → box with (s,u)
   - GLI["fctopology3", {1,1,1,1}] → crossed box with (t,u)
7. Full Born-virtual also includes vertex + self-energy + VP (not just box)

### Priority 2: MUnit test porting (~370 remaining)

High-value targets:
- feynfeld-8qe/n01/iaz/37v/32j: DiracTrick batches (~73 tests total)
- feynfeld-4mm: SUNSimplify (78 tests)
- feynfeld-36h: SUNTrace (24 tests)

### Priority 3 (backlog): Diagram generation epic (feynfeld-8c4)

Full FeynArts-class diagram generation in pure Julia. Four phases:
1. nauty/Traces for canonical graph labeling + topology enumeration
2. CSP-based field insertion with fermion flow + charge conservation
3. N-point amplitude mapping (generalize build_amplitude)
4. N-body phase space (RAMBO + VEGAS + Catani-Seymour dipoles)

~850 LOC total. Not in scope until NLO validation complete.

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
- **v2 source:** `src/v2/` (44 files, ~4,800 LOC)
- **v2 tests:** `test/v2/` (21 files + munit/) — 406 tests, ALL PASS
- **v1:** FROZEN. Do not extend.

### What Feynfeld can compute (textbook order)

**Tree-level QED:**
- e⁺e⁻ → μ⁺μ⁻ (massless, pipeline, P&S validated)
- Bhabha e⁺e⁻ → e⁺e⁻ (massless, 2 channels + interference, FeynCalc validated)
- Compton eγ → eγ (massive electron, P&S validated)

**Tree-level QCD:**
- qq̄ → gg (massless, 3 channels, SU(3) colour, FeynCalc validated)

**Tree-level EW:**
- e⁺e⁻ → W⁺W⁻ (massive W, 4 channels, chiral γ⁵, gauge cancellation, Grozin validated)

**1-loop QED (standalone PaVe):**
- Electron self-energy Σ(p²) via B₀, B₁, A₀ (Denner validated)
- Vertex correction F₂(0) = α/(2π) via C₀, C₁, C₂ (Schwinger validated)
- Vacuum polarization Π̂(q²), running α(q²) (PDG validated)
- Schwinger correction δσ/σ (analytical)
- Improved Born σ with VP + Schwinger

**1-loop QED (pipeline, NEW in Session 16):**
- Box e⁺e⁻ → μ⁺μ⁻: channels → amplitude → trace → TID → ComplexF64
- Direct box D₀(0,0,0,0,s,t) + crossed D₀(0,0,0,0,s,u)
- All 33 interference terms evaluate (rank 0/1/2 via PV + COLLIER)
- Awaiting FeynCalc cross-validation (feynfeld-73g)

**Cannot yet do:** φ⁴, 2→3+, full NLO pipeline (real emission, counter-terms),
2-loop, automatic diagram generation for arbitrary processes.

### Files new/modified in Session 16

| File | LOC | What |
|------|-----|------|
| `src/v2/loop_channels.jl` | 79 | NEW: LoopChannel + box_channels() |
| `src/v2/loop_amplitude.jl` | 115 | NEW: build_loop_box_amplitude() + BoxDenominators |
| `src/v2/loop_interference.jl` | 63 | NEW: spin_sum_tree_loop_interference() |
| `src/v2/tid.jl` | ~170 | NEW: evaluate_box_integral() + full rank-0/1/2 TID |
| `src/v2/pave_eval.jl` | +15 | MODIFIED: COLLIER B₀ fallback for B₀(0,0,0) |
| `src/v2/FeynfeldX.jl` | ~135 | MODIFIED: includes + exports |
| `test/v2/test_box_ee_mumu.jl` | 153 | NEW: 136 tests for Phases A-D |
| `test/v2/test_running_alpha.jl` | ~2 | MODIFIED: Δα sign fix |

### Box diagram momentum routing

Direct box e⁻(p₁)+e⁺(p₂)→μ⁻(k₁)+μ⁺(k₂):
```
D₀ = q²           (photon₁, V₂→V₁)
D₁ = (q+p₁)²      (electron, V₁→V₃)
D₂ = (q+p₁+p₂)²   (photon₂, V₃→V₄)
D₃ = (q+k₁)²       (muon, V₄→V₂)
PaVe: D₀(0,0,0,0,s,t, 0,0,0,0)
```
Crossed box: D₃ = (q+k₂)², PaVe: D₀(0,0,0,0,s,u, 0,0,0,0)

Accumulated momenta: K₁ = p₁, K₂ = p₁+p₂ (MomentumSum), K₃ = k₁ or k₂

### TID architecture

```
trace (AlgSum with SP(q,...)) → classify by rank → dispatch:
  rank 0: coeff × D₀(COLLIER)
  rank 1: Σᵢ (Kᵢ·pₐ) Dᵢ(PV reduction via d_tensor.jl)
  rank 2 SP(q,q): C₀⁽⁰⁾ + m₀²D₀
  rank 2 SP(q,pₐ)×SP(q,pᵦ): decompose pₐ in K basis (Gram inversion),
    then PV cancel → sub-triangle rank-1 (C₁/C₂ via PV, B₀ via COLLIER fallback)
```

### COLLIER interface notes

- Library: `refs/COLLIER/COLLIER-1.2.8/libcollier.so`
- Init: `ccall((:__collier_init_MOD_init_cll, lib), ...)` (module function, mangled name)
- D₀: `ccall((:d0_coli_, lib), ComplexF64, ...)` (plain COLI function)
- C₀: `ccall((:c0_coli_, lib), ComplexF64, ...)` (plain COLI function)
- B₀: `ccall((:b0_coli_, lib), ComplexF64, ...)` (plain COLI function, NEW in Session 16)
- D-tensor via module: `__collier_coefs_MOD_d_main_cll` exists but returns all zeros
  via ccall (assumed-shape array needs Fortran descriptors, not bare Ptr). Use PV instead.
- C-tensor via module: similar issue. Solved by using PV reduction + COLLIER B₀ fallback.

---

## QUICK COMMANDS

```bash
# Run single test (fast)
julia --project=. test/v2/test_box_ee_mumu.jl   # 5s, Spiral 10
julia --project=. test/v2/test_running_alpha.jl  # 4s, running α

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl            # 406 tests, ALL PASS

# Evaluate box integral at a kinematic point
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
# ... setup model, tree, loop amplitudes ...
result = evaluate_box_integral(numerator, sp_vals, denoms, s, t)
'

# Beads
bd ready              # available work
bd show feynfeld-73g  # NLO validation (NEXT)
bd show feynfeld-8c4  # diagram generation epic (BACKLOG)
bd stats              # project health

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
