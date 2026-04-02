# HANDOFF — 2026-04-02 (Session 15: eeZ γ5 fix, WWZ coupling sign, exact Grozin match)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (405+ tests, ~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 15 ACCOMPLISHMENTS

### 1. EXACT Grozin match: ee→WW at ALL energies (feynfeld-pyi + feynfeld-r3u CLOSED)

**The #1 priority from Session 14 is DONE.** Full e+e-→W+W- cross-section now
matches the Grozin/Denner analytical formula to ratio = 1.0000 at all LEP2 energies.

| √s (GeV) | σ_pipeline (pb) | σ_Grozin (pb) | Ratio |
|-----------|-----------------|---------------|-------|
| 170 | 14.252 | 14.252 | 1.0 |
| 200 | 18.01 | 18.01 | 1.0 |
| 300 | 12.556 | 12.556 | 1.0 |
| 500 | 6.658 | 6.658 | 1.0 |

**Two bugs found and fixed simultaneously:**

**Bug 1: eeZ vertex γ5 ordering** (`rules.jl:27-31`)
The vertex had `(gV - gA γ5)γ^μ` (γ5 BEFORE γ^μ) instead of the standard
`γ^μ(gV - gA γ5)` (γ5 AFTER γ^μ). Since `γ5 γ^μ = -γ^μ γ5`, the wrong
ordering gave `γ^μ(gV + gA γ5)`, effectively swapping L↔R chirality.
This flipped the sign of gA in the sZ×tν cross-term.

Fix: `DiracChain([GA5(), γ^μ])` → `DiracChain([γ^μ, GA5()])`.
Ref: Denner1993 Eq. (A.13), PDG2024 Table 10.3.

**Bug 2: WWZ coupling phase** (`rules.jl:48-53`, `test_ee_ww_grozin.jl:39`)
The VVV coupling C (Denner Eq. A.7) is C=+1 for γWW but C=-cW/sW for ZWW.
The SIGN of C creates a relative phase between M_sγ and M_sZ amplitudes.
The coupling weight `const_c_sZ` was positive (magnitude only), missing the
negative sign from C_ZWW. This caused the sZ×tν cross-term to have the
wrong sign relative to the sγ×tν cross-term.

Fix: Added `gauge_coupling_phase(Val(:g_WWZ)) = -1//1` dispatch function.
Applied in test: `const_c_sZ = gauge_coupling_phase(Val(:g_WWZ)) * e²/(2sin²θ_W)`.
Ref: Denner1993 Eq. (A.7): "ZW⁺W⁻: C = -c/s".

**Why these bugs compensated each other:** With the old vertex (Bug 1, gA→-gA
in cross-terms) AND the positive coupling weight (Bug 2, missing -1 on sZ),
the product of two wrong signs gave a partially correct cross-term. Near
threshold (170 GeV) the error was 2%; at 500 GeV it grew to 55%.

### 2. Three additional bug fixes

**`_scalar` helper** (test file): `first(r.terms)` picked an arbitrary Dict
entry, giving phantom values when Eps terms were present. Fixed to use
`FactorKey()` lookup for the true scalar coefficient.

**`fermion_contract` break condition** (test file): Missing Eps acceptance
in the convergence check (ran all 10 iterations instead of breaking early).
Fixed to match `gauge_contract`'s condition.

**`propagator_num(::Boson)`**: Identified as dead code — defined but never
called. Not fixed (no functional impact), but documented for future cleanup.

### 3. Beads

- **Closed:** feynfeld-pyi (eeZ γ5 ordering bug)
- **Closed:** feynfeld-r3u (deep investigation of cross-term sign)

---

## KNOWN ISSUES AND BLOCKERS

### P1 (pre-existing): Δα imaginary part — 1 flaky test
Intermittent. Not related to Session 15 changes.

### P3: `propagator_num(::Boson)` is dead code
Defined in rules.jl but never called. The boson propagator numerator -g_{μν}
is not applied in `_build_gauge_exchange` or `_build_boson_exchange`. For
now this is harmless because the (-1)² cancels in diagonals and the relative
sign is handled by `gauge_coupling_phase`. Future: consider activating it.

### P3: `PropagatorRule` exported but undefined
Pre-existing: FeynfeldX.jl exports `PropagatorRule` but the struct doesn't exist.

---

## WHAT TO DO NEXT

### Priority 1: MUnit test porting (~370 remaining)

MUnit progress: 60/~430 done across 5 function files + 2 batch files.
High-value next targets:
- feynfeld-rcd: Contract D-dim section (20 tests, all portable)
- feynfeld-4mm: SUNSimplify (78 tests)
- feynfeld-36h: SUNTrace (24 tests)

### Priority 2: Spiral 10 continuation
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
- **v2 source:** `src/v2/` (39 files, ~4,200 LOC)
- **v2 tests:** `test/v2/` (20 files + munit/) — 405+ tests
- **v1:** FROZEN. Do not extend.

### New/modified files in Session 15

| File | LOC | What |
|------|-----|------|
| `src/v2/rules.jl` | 116 | eeZ vertex fix + gauge_coupling_phase dispatch |
| `src/v2/FeynfeldX.jl` | ~140 | Export gauge_coupling_phase |
| `test/v2/test_ee_ww_grozin.jl` | 289 | Coupling sign fix, _scalar fix, tolerance tightened |

### Coupling constant conventions (UPDATED from Session 14)

| Channel | Pipeline output | Coupling weight | Propagator |
|---------|----------------|-----------------|------------|
| s-γ | DiracExpr[γ^ρ] + AlgSum[V_ρμν] | +e² = +4πα | 1/s |
| s-Z | DiracExpr[(gV−gAγ5)γ^ρ] + AlgSum[V_ρμν] | **−e²/(2sin²θ_W)** | 1/(s−M_Z²) |
| t-ν | DiracChain[γ^μ P_L q̸ γ^ν P_L] | +e²/(2sin²θ_W) | 1/t |

Note: **c_sZ is NEGATIVE** (from gauge_coupling_phase). This encodes the relative
phase between γWW (C=+1) and ZWW (C=−cW/sW) triple gauge couplings (Denner A.7).
The sign cancels in diagonals (w_Z² > 0) but matters for interference terms.

---

## QUICK COMMANDS

```bash
# Run single test (fast)
julia --project=. test/v2/test_vertical.jl       # 5s, pipeline
julia --project=. test/v2/test_ee_ww_grozin.jl   # 12s, Stage C validation (ratio = 1.0)
julia --project=. test/v2/test_pipeline.jl        # 8s, all processes

# Run MUnit tests
for f in test/v2/munit/test_*.jl; do julia --project=. "$f"; done

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Beads
bd ready              # available work
bd stats              # project health

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
