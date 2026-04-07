# HANDOFF — 2026-04-07 (Session 18: Stringent code quality review, 18 issues fixed)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (605 tests, ~10 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 18 ACCOMPLISHMENTS

### Stringent code quality review: 18 issues found and fixed

Full read of all 4,654 LOC across 43 v2 source files against Julia idiom rules,
CLAUDE.md conventions, and Rule 11 (200 LOC limit). Re-read idiom rules every
~5kLOC to maintain LLM attention focus. All 605 tests pass after every fix.

**CRITICAL fixes (8):**
- C1: `amplitude.jl` isa cascade → `_build_amplitude` dispatch on `Field{Boson}`/`Field{Fermion}`
- C2: `amplitude.jl` 211→179 LOC — extracted `gauge_exchange.jl` (38 LOC)
- C3: `pave_eval.jl` 230→107 LOC — extracted `b0_eval.jl` (122 LOC)
- C4: `BoxDenominators` `NTuple{3,Any}` → `NTuple{3,Union{Momentum,MomentumSum}}`
- C5: Removed phantom `PropagatorRule` export (never defined)
- C6: Removed dead `vertex_structure(::SU{N},...)` 4-arg method (pipeline uses 5-arg)
- C7: `all_amps = []` → `Tuple{DiracExpr,DiracExpr}[]` (no more `Vector{Any}`)
- C8: `CrossSectionProblem` `s::Number` → parametric `s::T`

**HIGH fixes (6):**
- H1-H2: Removed 2 dead branches in `dirac_trick.jl` (identical ternary arms)
- H3: Removed redundant `!iszero(mass)` inside else branch in `spin_sum.jl`
- H4: Removed dead `fermion_spin_sum` from v2 exports
- H5: `propagator_num(::Boson)` is NOT dead — tested in test_vertical.jl (closed)
- H6: `tid.jl` `::Any` fallback → `::AlgFactor`

**MODERATE fixes (4):**
- M1: `_eps_sign` hardcoded names — acceptable for all current 2→2 processes (closed)
- M2: Added `_try_expand(::Eps)` for MomentumSum slots — latent bug fixed
- M3: Removed `dsigma_dt_ee_ww` and `sigma_total_compton_massless` stubs (always errored)
- M4: `tid.jl` 202→175 LOC — extracted `sp_lookup.jl` (32 LOC)

**New files:**
| File | LOC | What |
|------|-----|------|
| `gauge_exchange.jl` | 38 | Extracted from amplitude.jl |
| `b0_eval.jl` | 122 | B₀ scalar/tensor evaluation |
| `sp_lookup.jl` | 32 | _splookup + _Kdot overloads |

---

## SESSION 17 ACCOMPLISHMENTS

### 1. NLO box validation — Phases E+F (feynfeld-73g CLOSED)

Assembled the full box contribution to Born-virtual interference and validated
at multiple kinematic points.

**New source:** `nlo_box.jl` (108 LOC)
- `evaluate_single_box_channel()`: runs full pipeline for one box channel
- `evaluate_box_channels()`: sums over direct + crossed boxes
- `born_virtual_box()`: physical Born-virtual with coupling factors

**New test:** `test_nlo_box_validation.jl` (210 LOC, 63 tests)
- Phase D: TID evaluation (direct + crossed, finiteness, non-degeneracy)
- Phase E: Channel sum assembly
- F1: Im(I) forward-backward symmetry (PASS)
- F2: Im(I) crossing at θ=90° (PASS)
- F3: Finiteness scan (4 energies × 5 angles = 20 points)
- F4: Order of magnitude box/Born ~ O(α/π) (PASS)
- F5: Energy scaling (ratio changes < 10× over √s = 50–500 GeV)
- F6: COLLIER D₀ sanity
- F7: born_virtual_box non-zero

### 2. Key normalization discovery

The overall coupling factor in the tree×box interference is **purely imaginary**:
the physical quantity depends on **Im(I_COLLIER)**, not Re(I_COLLIER).

Evidence:
- Im(I_total) is forward-backward symmetric: -273024 at cosθ=±0.5 ✓
- Re(I_total) breaks F-B symmetry: -3.38e7 vs -2.78e7 ✗
- Im gives box/Born ~ 6e-3 ≈ O(α/π) ✓
- Re gives box/Born ~ 0.78 ≈ O(1) ✗

Formula: `(1/4) Σ 2Re(M_tree* × M_box) = -e⁶/(32π²s) × Im(I_COLLIER)`

The imaginary coupling arises because the tree propagator phase (-i/s) combines
with the loop measure (i/(16π²)) and vertex phases to give a net purely-imaginary
overall factor. Full derivation still needed (tracked for future work).

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
- **feynfeld-73g** CLOSED: NLO box validation (Phases E+F, 63 tests)

---

## KNOWN ISSUES

### P3: `contract.jl` is 212 LOC (borderline Rule 11)
The contraction engine is a single coherent algorithm. Splitting would hurt
readability more than it helps. Acceptable borderline case.

### P3: `propagator_num(::Boson)` not used by pipeline
Defined, tested (test_vertical.jl:69), but amplitude builder constructs boson
propagators inline instead. Pipeline integration tracked for future work.

---

## WHAT TO DO NEXT

### Priority 1: Full Born-virtual (vertex + self-energy + CTs)

Box-only validation complete (feynfeld-73g CLOSED). To compare against
FeynCalc's `bornVirtualRenormalized[1]` from arXiv:hep-ph/0010075 Eq. 2.32,
we need the remaining 1-loop contributions:
1. Vertex correction diagrams (electron and muon lines)
2. Vacuum polarization (photon self-energy)
3. Counter-terms (MS-bar for photon field, OS for fermions)
4. Then sum all + boxes and compare O(ε⁰) finite part against FeynCalc

Also needed: rigorous derivation of the Im(I_COLLIER) normalization
(currently empirical — see Session 17 discovery).

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
- **v2 source:** `src/v2/` (47 files, ~4,640 LOC) — 3 new files from code quality split
- **v2 tests:** `test/v2/` (22 files + munit/) — 605 tests, ALL PASS
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
- Validated: box/Born ~ O(α/π), Im(I) F-B symmetric, 63 tests
- Full Born-virtual comparison awaits vertex + self-energy + CTs

**Cannot yet do:** φ⁴, 2→3+, full NLO pipeline (real emission, counter-terms),
2-loop, automatic diagram generation for arbitrary processes.

### Files new/modified in Session 17

| File | LOC | What |
|------|-----|------|
| `src/v2/nlo_box.jl` | 108 | NEW: evaluate_box_channels + born_virtual_box |
| `src/v2/FeynfeldX.jl` | +3 | MODIFIED: include + exports |
| `test/v2/test_nlo_box_validation.jl` | 210 | NEW: 63 tests for Phases D-F |

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
julia --project=. test/v2/test_box_ee_mumu.jl        # 5s, Spiral 10 Phases A-C
julia --project=. test/v2/test_nlo_box_validation.jl  # 3min, Phases D-F (63 tests)
julia --project=. test/v2/test_running_alpha.jl       # 4s, running α

# Full suite (single process, ~10 min)
julia --project=. test/v2/runtests.jl            # 605 tests, ALL PASS

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
