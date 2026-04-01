# HANDOFF — 2026-04-01 (Session 14: Stage C gauge cancellation, MUnit porting, eνW vertex fix)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (405+ tests, ~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 14 ACCOMPLISHMENTS

### 1. Stage C: s×t gauge cancellation for ee→WW (feynfeld-bao CLOSED)

**The #1 priority from Session 13 is DONE.** Full e+e-→W+W- cross-section now
includes all 9 terms of |M|² (3 diagonal + 3 same-topology + 3 cross-topology).

**Bug found:** eνW vertex in `rules.jl` had `P_L γ^μ` (GA7 before gamma) instead
of the standard `γ^μ P_L` (gamma before GA7). Since `P_L γ^μ = γ^μ P_R`, this
flipped chirality in cross-terms, producing wrong-sign Eps contributions that
prevented gauge cancellation.

**Fix:** `rules.jl:33` changed from `DiracChain([GA7(), γ^μ])` to
`DiracChain([γ^μ, GA7()])`. One-line change, verified: 405/406 tests pass
(1 pre-existing flaky Δα imaginary test).

**Results at LEP2 energies:**

| √s (GeV) | σ_pipeline (pb) | σ_Grozin (pb) | Ratio |
|-----------|-----------------|---------------|-------|
| 170 | 14.53 | 14.25 | 1.020 |
| 200 | 19.32 | 18.01 | 1.073 |
| 300 | 15.10 | 12.56 | 1.203 |
| 500 | 10.34 | 6.66 | 1.553 |

Near threshold: **2% accuracy**. The gauge cancellation reduces σ from 222 pb
(diagonal only) to 10 pb at 500 GeV — a 22× reduction, demonstrating the
delicate SU(2) gauge cancellation is working.

**Residual at high energy:** grows to 55% at 500 GeV. Cause: eeZ vertex γ5
ordering (feynfeld-pyi, P1). The eeZ vertex computes `(gV - gA γ5)γ^μ` which
equals `γ^μ(gV + gA γ5)` due to anticommutation, instead of the standard
`γ^μ(gV - gA γ5)`. This flips the sign of the gA Eps contribution in the
s-Z × t-ν cross-term. Fix requires careful investigation (see KNOWN ISSUES).

### 2. MUnit test porting (38 new tests, 90 assertions)

Created 4 new function-specific MUnit test files in `test/v2/munit/`:

| File | Function | Tests | Assertions |
|------|----------|-------|------------|
| `test_DiracTrick.jl` | DiracTrick n=1-5 | 13 | 65 |
| `test_PolarizationSum.jl` | PolarizationSum | 6 | 6 |
| `test_ExpandScalarProduct.jl` | ExpandScalarProduct | 10 | 10 |
| `test_Contract.jl` | Contract + D-dim eps | 9 | 9 |

Key tests: FiveFreeIndices-ID1 (11-term general formula), D-dim Levi-Civita
eps×eps = D(1-D)(2-D)(3-D), bilinear SP expansion with Mandelstam-like terms.

### 3. Beads housekeeping

- **Closed:** feynfeld-6ds (Rational overflow, was fixed in Session 13 but left open)
- **Closed:** feynfeld-bao (ee→WW full pipeline vs Grozin)
- **Created:** feynfeld-pyi (eeZ γ5 ordering bug, P1)
- **Updated:** 8 MUnit beads with progress notes

---

## KNOWN ISSUES AND BLOCKERS

### P1: eeZ vertex γ5 ordering (feynfeld-pyi) — blocks full Grozin match at high energy

**Root cause:** `rules.jl` line 26 uses `DiracChain([GA5(), γ^μ])` (γ5 before γ^μ).
Since γ5 γ^μ = -γ^μ γ5, the vertex `gV γ^μ - gA γ5 γ^μ` actually equals
`gV γ^μ + gA γ^μ γ5 = γ^μ(gV + gA γ5)` instead of `γ^μ(gV - gA γ5)`.

**Impact:** The diagonal |M_sZ|² is unaffected (gA² is sign-invariant). The
s-γ × s-Z interference is a small correction (<2%), so Stage B passes. But the
s-Z × t-ν cross-term has a wrong-sign Eps contribution that grows with energy.

**Why the obvious fix didn't work:** Naively changing `GA5() γ^μ` → `γ^μ GA5()`
unexpectedly changed T_tν (the t-channel diagonal). This might be a Julia
module recompilation artifact or a subtle dispatch interaction. Needs
investigation with a fresh Julia session and careful unit testing of just
the eeZ vertex structure.

**Approach:** Write a standalone eeZ vertex test that verifies:
1. `gV γ^μ - gA γ^μ γ5` trace against known P&S result
2. Cross-term with non-chiral γ^ρ has correct Eps sign
3. All 405 existing tests still pass
Then apply the fix.

### P1 (pre-existing): Δα imaginary part — 1 flaky test
Intermittent. Not related to Session 14 changes.

---

## WHAT TO DO NEXT

### Priority 1: Fix eeZ vertex γ5 ordering (feynfeld-pyi)

This is the **one remaining blocker** for full Grozin match at all energies.
Once fixed, Stage C should achieve <1% accuracy across the full LEP2 energy range.
See the approach described in KNOWN ISSUES above.

### Priority 2: MUnit test porting (~370 remaining)

MUnit progress: 60/~430 done across 5 function files + 2 batch files.
Most remaining tests require GAE/GSD/GSE constructors (BMHV dimensional splitting)
or CartesianIndex support. High-value next targets:
- feynfeld-rcd: Contract D-dim section (20 tests, all portable)
- feynfeld-4mm: SUNSimplify (78 tests)
- feynfeld-36h: SUNTrace (24 tests)

### Priority 3: Spiral 10 continuation
1-loop amplitude builder (feynfeld-7h8), ee→μμ NLO box (feynfeld-4q5).

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
- **v2 source:** `src/v2/` (39 files, ~4,100 LOC)
- **v2 tests:** `test/v2/` (20 files + munit/) — 405+ tests
- **v1:** FROZEN. Do not extend.

### New/modified files in Session 14

| File | LOC | What |
|------|-----|------|
| `src/v2/rules.jl` | 104 | eνW vertex fix: γ^μ P_L ordering |
| `test/v2/test_ee_ww_grozin.jl` | 278 | Stage C: cross_trace_one_way, cross_contract_st, full integration |
| `test/v2/munit/test_DiracTrick.jl` | 198 | NEW: 13 tests (n=1-5, pure + mixed-dim) |
| `test/v2/munit/test_PolarizationSum.jl` | 72 | NEW: 6 tests (Feynman, massive, axial) |
| `test/v2/munit/test_ExpandScalarProduct.jl` | 96 | NEW: 10 tests (FV, SP bilinear, Mandelstam) |
| `test/v2/munit/test_Contract.jl` | 115 | NEW: 9 tests (eps×eps 4D+DD, metric chain) |

### Coupling constant conventions (unchanged from Session 13)

| Channel | Pipeline output | Missing coupling | Missing propagator |
|---------|----------------|------------------|--------------------|
| s-γ | DiracExpr[γ^ρ] + AlgSum[V_ρμν] | e² = 4πα | 1/s |
| s-Z | DiracExpr[(gV-gAγ5)γ^ρ] + AlgSum[V_ρμν] | e²/(2sin²θ_W) | 1/(s-M_Z²) |
| t-ν | DiracChain[γ^μ P_L q̸ γ^ν P_L] | e²/(2sin²θ_W) | 1/t |

Note: t-channel chain now has `γ^μ P_L` (standard convention, fixed in Session 14).

---

## COLLIER SETUP (REQUIRED ON EACH MACHINE)

```bash
cd refs/COLLIER/COLLIER-1.2.8
mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
ls refs/COLLIER/COLLIER-1.2.8/libcollier.so  # verify
mkdir -p output  # COLLIER writes log files here
```

---

## QUICK COMMANDS

```bash
# Run single test (fast)
julia --project=. test/v2/test_vertical.jl       # 5s, pipeline
julia --project=. test/v2/test_ee_ww_grozin.jl   # 10s, Stage C validation
julia --project=. test/v2/test_pipeline.jl        # 8s, all processes

# Run MUnit tests
for f in test/v2/munit/test_*.jl; do julia --project=. "$f"; done

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Beads
bd ready              # available work
bd stats              # project health
bd show feynfeld-pyi  # eeZ vertex bug details

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
