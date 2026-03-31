# HANDOFF — 2026-03-31 (Session 13: Rational overflow fix, Eps contraction, Phase B)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (405+ tests, ~5 min)
4. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 13 ACCOMPLISHMENTS

### 1. Fixed Rational{Int} overflow (feynfeld-6ds)
**Root cause:** `Mandelstam(s, cosθ)` called `Rational{Int}(cosθ)` when cosθ was Float64
from quadgk. Julia converts Float64 → Rational with enormous denominators → overflow in
subsequent `mul_coeff` calls.

**Fix:** Parametric `Mandelstam{T<:Real}` + Float64 evaluation path that bypasses SPContext.

**Files changed:**
- `src/v2/cross_section.jl`: Mandelstam{T}, evaluate_numeric, sp_values_2to2
- `src/v2/eps_evaluate.jl` (NEW): dispatch-based factor evaluation, Eps Gram determinant
- `src/v2/FeynfeldX.jl`: includes + exports

**Design:** Symbolic algebra stays purely Rational{Int}. Float64 enters only at the
evaluation boundary. The pipeline produces a symbolic AlgSum, then `evaluate_numeric`
walks it with Float64 arithmetic. Mirrors FeynCalc (symbolic until `N[]`) and
FormCalc (symbolic → compiled Fortran).

### 2. Eps contraction engine fix
**Root cause:** `contract.jl` had NO dispatch methods for Eps × MetricTensor or
Eps × FourVector. Chiral Z coupling produces γ5 traces → Eps tensors that couldn't
be contracted with the triple gauge vertex's metric/four-vector factors.

**Fix:** Added `_do_contraction(::Eps, ::MetricTensor, ...)` and
`_do_contraction(::Eps, ::FourVector, ...)` to contract.jl, plus antisymmetric
vanishing detection (`_eps_replace_slot` returns `:vanishes` when two slots become identical).

**Citations:**
- Eps × Metric/FV: FeynCalc PairContract.m lines 169-170; EpsEvaluate.m lines 107-111
- Gram determinant: MertigBohmDenner1991 Eq. (2.21)
- Antisymmetric vanishing: FeynCalc EpsEvaluate.m `Signature[{x}] === 0`

### 3. Grozin formula validated against Denner (1993)
`sigma_ee_ww()` matches Denner Tab. 11.4 Born(G_F) to 0.01% at all LEP2 energies
when using Denner's parameters (M_W=80.23, α from G_F).

### 4. Multi-channel ee→WW validation (Stages A+B)
`test/v2/test_ee_ww_grozin.jl` validates all 3 diagonal |M_i|² (s-γ, s-Z, t-ν)
plus s-γ × s-Z interference. Results are physically correct:
- Diagonal sum exceeds Grozin total (ratio 1.4-33, growing with energy)
- γ-Z interference is destructive (reduces σ by ~3-5%)
- Remaining gap is the s×t gauge cancellation (Stage C)

---

## KNOWN ISSUES AND BLOCKERS

### P1: Stage C — s×t cross-topology interference (blocks full Grozin match)
The s-channel (gauge exchange) and t-channel (fermion exchange) have different
topologies. Their interference requires a cross-trace combining:
- s-channel: single γ^ρ vertex → rho_s connects to triple gauge vertex
- t-channel: P_L γ^μ q-slash γ^ν P_L → mu_k1, mu_k2 connect to W polarizations

After the cross-trace, the ρ_s index needs contraction with the triple gauge vertex
while μ indices contract with polarization sums. This is a SINGLE fermion line
(not reconnected), so it's a standard Dirac trace, not a Bhabha-style interference.

**Approach:** Build a DiracExpr that sums s-channel and t-channel fermion structures
(different gamma sequences but same spinors). The s-channel terms have ρ_s index,
the t-channel terms have μ_k1/μ_k2. After _single_line_trace, cross-terms have mixed
indices. Contract the s-channel index with triple_gauge_vertex and polarization sums.

**This is the dominant gauge cancellation** — without it, σ grows as s²/M_W⁴ instead
of log(s)/s (Denner Eqs. 11.16 vs 11.17).

### P2: Eps contraction engine — multi-pass requirement (feynfeld-qu1)
Large chiral expressions (s-Z channel: 1024 initial terms) require up to 10
contraction passes to fully reduce. The worklist algorithm in `_contract_factors`
finds only one contractible pair per pass. This is correct but slow. A future
optimization could batch-contract all pairs found in one scan.

### P1 (pre-existing): Δα imaginary part (timelike) — 1 flaky test
Running α test in `test/v2/test_running_alpha.jl`, "Δα imaginary part (timelike)"
section. Fails intermittently. Not related to Session 13 changes.

---

## WHAT TO DO NEXT

### Priority 1: Stage C — s×t interference for full Grozin match

**This is the #1 task.** Everything is in place:
- All 3 channels build through pipeline ✓
- Diagonal |M_i|² evaluate numerically ✓
- s-γ × s-Z interference works ✓
- Float64 integration over cosθ works ✓
- Coupling constants and propagators verified ✓

**What's needed:** Add the s×t cross-terms to `test/v2/test_ee_ww_grozin.jl`.

**Concrete approach (validated by research):**

All three channels share ONE fermion line (incoming e+e-). The s-channel fermion
chain has elements `[v̄, γ^{ρ_s}, u]` (or chiral variant), while the t-channel has
`[v̄, P_L γ^{μ_k1} GS(q) γ^{μ_k2} P_L, u]`. These have the SAME spinors.

Create a combined DiracExpr with terms from all 3 channels (1 from s-γ, 2 from s-Z,
1 from t-ν = 4 terms total). Call `_single_line_trace` to get ALL 16 (i,j) pairs.

After tracing, cross-terms have MIXED indices:
- s×s terms: only ρ_s, ρ_s_ → contract with vtx × vtx_conj × P1 × P2
- t×t terms: only μ_k1, μ_k2, μ_k1_, μ_k2_ → contract with P1 × P2
- s×t terms: ρ_s + μ_k1_, μ_k2_ (or vice versa) → contract with vtx on s-side, P on both

The challenge: after tracing, the mixed-index terms need DIFFERENT contraction
treatments for the s-side and t-side indices. This may require:
1. Separate extraction of cross-terms (T_combined - T_s_only - T_t_only)
2. Contracting with the appropriate vertex + polarization structure
3. Careful index management

**Alternative simpler approach:** Compute the s×t cross-trace DIRECTLY by building
the trace manually: `Tr[completeness × Γ_s × completeness × Γ̄_t]` where Γ_s and Γ_t
are the gamma sequences from each channel. Then contract with the appropriate external
structures. This avoids the combined-DiracExpr index-mixing problem.

**Validation:** The result (diagonal + all interference) should match `sigma_ee_ww()`
to high precision at all energies. The ratio should be ~1.0 ± numerical integration
tolerance (~0.1%).

### Priority 2: MUnit test porting (~386 remaining)
See Session 12 HANDOFF for the full bead list. 22/408 done (DiracTrace).

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
9. **NEVER modify TensorGR.jl or core algebra files without explicit permission.**
   Session 13 learned this the hard way — got caught modifying contract.jl without
   following the tiered workflow. DO NOT repeat.

---

## PROJECT STATE

### Branch and code location
- **Branch:** `master`
- **v2 source:** `src/v2/` (38 files, ~4,100 LOC)
- **v2 tests:** `test/v2/` (20 files + munit/) — 405+ tests
- **v1:** FROZEN. Do not extend.

### New/modified files in Session 13

| File | LOC | What |
|------|-----|------|
| `src/v2/cross_section.jl` | 156 | Mandelstam{T}, evaluate_numeric, sp_values_2to2 |
| `src/v2/eps_evaluate.jl` | 65 | NEW: _eval_factor dispatch, _evaluate_eps Gram det |
| `src/v2/contract.jl` | 212 | Eps×Metric, Eps×FV contraction + vanishing |
| `src/v2/FeynfeldX.jl` | 121 | +2 lines: include + export |
| `test/v2/test_ee_ww_grozin.jl` | 163 | NEW: Stages A+B validation |

### Coupling constant conventions (CRITICAL for Stage C)

Pipeline `build_amplitude` returns Lorentz/Dirac structure ONLY. Missing factors:

| Channel | Pipeline output | Missing coupling | Missing propagator |
|---------|----------------|------------------|--------------------|
| s-γ | DiracExpr[γ^ρ] + AlgSum[V_ρμν] | e² = 4πα | 1/s |
| s-Z | DiracExpr[(gV-gAγ5)γ^ρ] + AlgSum[V_ρμν] | e²/(2sin²θ_W) | 1/(s-M_Z²) |
| t-ν | DiracChain[P_L γ^μ q̸ γ^ν P_L] | e²/(2sin²θ_W) | 1/t |

Ref: Denner1993, Eq. (11.9); PDG2024, Table 10.3.
g_V, g_A are IN the s-Z vertex_structure. P_L is IN the t-ν chain.
All other couplings must be applied externally.

In M_W units: M_Z² = 1/cos²θ_W ≈ 1.288. All masses in M_W² = 1.

### Trace index conventions

| Channel | Trace indices (forward) | Trace indices (conjugate) |
|---------|------------------------|--------------------------|
| s-γ | rho_s | rho_s_ |
| s-Z | rho_s | rho_s_ |
| t-ν | mu_k1, mu_k2 | mu_k1_, mu_k2_ |

After _single_line_trace, forward indices connect to vtx/P1/P2 (original),
conjugate indices connect to vtx_conj/P1/P2 (conjugate). For s×t cross-terms,
mixed indices appear (rho_s from s-side, mu_k1_ from t-side conjugate).

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
julia --project=. test/v2/test_ee_ww_grozin.jl   # 9s, multi-channel validation
julia --project=. test/v2/test_pipeline.jl        # 8s, all processes

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Beads
bd ready              # available work
bd stats              # project health
bd show feynfeld-qu1  # Eps contraction bug details

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
