# HANDOFF — 2026-03-29 (End of Session 10, COLLIER integration)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Read `Feynfeld_PRD.md` — vision, endgame interaction, spiral plan, hardness scale
3. Read `src/v2/DESIGN.md` — type system, anti-patterns, Session 8 review findings
4. Run `bd ready` to see available work
5. **BUILD COLLIER** (see section below — required on each machine)
6. Run `julia --project=. test/v2/test_vertical.jl` to verify pipeline works
7. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## THE PIPELINE PRINCIPLE

**The pipeline IS the architecture. If a process bypasses it, the architecture is incomplete.**

Every process MUST flow through: Model → Rules → Channels → Amplitude → Algebra → Evaluate.
No hand-built amplitudes in test files. No standalone recipes that skip the pipeline.
Standalone formulas (schwinger.jl, ew_cross_section.jl, etc.) are REFERENCE IMPLEMENTATIONS
for cross-validation — not substitutes for the pipeline.

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **CITE EVERYTHING (tiered).**
3. **ACQUIRE PAPERS BEFORE CODE.** arXiv → TIB VPN → playwright-cli.
4. **JULIA IDIOMATIC.** Dispatch, not isa cascades. No Any.
5. **WORKFLOW (TIERED):** <5 LOC direct, <20 LOC 1+1, >20 LOC 3+1.
6. **REVIEW.** Rigorous reviewer after every core change.
7. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.
8. **LOC LIMIT ~200.** No source file exceeds ~200 lines.
9. **NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT OVERVIEW

### What is Feynfeld.jl?

Julia-native, agent-facing, full-stack physics computation suite. Lagrangian →
cross-section in one `using Feynfeld`. Replaces FeynRules + FeynArts + FeynCalc +
FormCalc + LoopTools. See PRD for the full vision.

### Branch and code location

- **Branch:** `master` (experimental/rebuild-v2 merged into master Session 10)
- **v2 source:** `src/v2/` (35 files, ~3,600 LOC)
- **v2 tests:** `test/v2/` (17 files, 329 tests)
- **v1:** `src/algebra/`, `src/integrals/` — FROZEN, will be deleted. Do NOT extend.

---

## COLLIER SETUP (REQUIRED ON EACH MACHINE)

COLLIER is the scalar loop integral library. `libcollier.so` is NOT checked in.
If missing, C₀ falls back to slow quadgk (no crash — graceful fallback).

```bash
# Download and build (one time per machine)
cd refs/COLLIER/COLLIER-1.2.8
mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
cd /path/to/Feynfeld.jl
ls refs/COLLIER/COLLIER-1.2.8/libcollier.so  # verify

# COLLIER writes log files here
mkdir -p output
```

**If COLLIER source not present:** Download from HepForge:
```bash
mkdir -p refs/COLLIER
curl -sL "https://www.hepforge.org/archive/collier/collier-1.2.8.tar.gz" | tar xz -C refs/COLLIER
```

---

## WHAT WAS DONE IN SESSION 10

### 1. Merged experimental/rebuild-v2 → master
Deleted experimental branch locally and on remote. All development on `master`.

### 2. COLLIER integration for C₀ (the main achievement)

**File:** `src/v2/c0_analytical.jl` (82 LOC)

Replaced broken Denner/Spence analytical C₀ with ccall to COLLIER (GPL v3).
Three evaluation paths:
- **C0p0** (all p²=0): analytical, pure logarithms
- **General**: COLLIER `C0_coli` ccall (~μs, handles spacelike/timelike/threshold)
- **Gram-degenerate** (Δ₂=0): quadgk fallback (~30s, rare edge case)

**Key results:** C₁/C₂ tensor reduction went from ~60s to 0.1s (600× speedup).
Spacelike C₀ went from ~30s to ~μs (100,000× speedup).

### 3. Lessons learned (CRITICAL — read before touching loop integrals)

- **Denner C0p3 formula (C0func.F) is a BACKUP, not the primary algorithm.**
  It has broken η corrections for Kallen < 0. The primary is FF (ffxc0).
  We wasted hours on this. DO NOT try to fix the Denner formula.

- **LoopTools `spence(0, x, 0)` = Li₂(x)**, NOT Li₂(1-x). Verified by reading
  auxCD.F: Li2series(z1) computes Li₂(1-z1), and spence(0,x,0) calls
  Li2series(1-x) = Li₂(x).

- **LoopTools.jl returns WRONG values** for Kallen < 0 unless you call
  `LoopTools.setwarndigits(100)` first. Default warndigits=9 causes FF
  result to be overridden by broken Denner backup.

- **Gram-degenerate scalar C₀** (Δ₂=0) cannot be computed by ANY Spence-based
  algorithm — the 1/√Δ₂ factor diverges. Neither COLLIER nor LoopTools
  handles it. quadgk is correct.

- **COLLIER argument convention:** Feynfeld (p10,p12,p20,m02,m12,m22) maps
  directly to COLLIER C0_coli args — same ordering, NO reorder.

- **COLLIER ccall pattern:**
  ```julia
  # Init (once, requires mkdir -p output):
  ccall((:__collier_init_MOD_init_cll, lib), Cvoid,
        (Ref{Int32}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Csize_t),
        Ref{Int32}(4), Ref{Int32}(4), "output", Ref{Int32}(0), 6)

  # C0 (plain Fortran function, symbol c0_coli_):
  ccall((:c0_coli_, lib), ComplexF64,
        (Ref{ComplexF64}, ..., Ref{ComplexF64}),
        ComplexF64(p10), ComplexF64(p12), ComplexF64(p20),
        ComplexF64(m02), ComplexF64(m12), ComplexF64(m22))
  ```

### 4. Papers acquired
- `refs/papers/vanOldenborgh1990_ZPhysC46.pdf` — FF algorithms paper

### 5. LoopTools.jl added as dependency
In Project.toml. Useful as test oracle (with setwarndigits fix).

---

## WHAT TO DO NEXT

### Priority 1: D₀ via COLLIER ccall (~30 min)

**THE OBVIOUS NEXT STEP.** Same pattern as C₀. COLLIER has `d0_coli_`:
```bash
nm -D refs/COLLIER/COLLIER-1.2.8/libcollier.so | grep d0_coli
```
D₀ takes 10 arguments: 6 momentum invariants + 4 masses. Check signature in
`refs/COLLIER/COLLIER-1.2.8/src/COLI/coli_d0.F`. Unblocks Spiral 10 (box diagrams).

| Issue | What | Effort |
|-------|------|--------|
| `feynfeld-62k.5` | D₀ via COLLIER ccall | ~30 min |

### Priority 2: Single-process test runner (~1 hour)

| Issue | What | Impact |
|-------|------|--------|
| `feynfeld-icg` | `test/v2/runtests.jl` — one Julia process | 12 min → 2 min |

### Priority 3: Cross-validate PaVe against COLLIER (~1 hour)

| Issue | What | Impact |
|-------|------|--------|
| `feynfeld-62k.9` | Systematic A₀/B₀/C₀ comparison | Trust in numerics |

### Priority 4: Pipeline completion

| Issue | What |
|-------|------|
| ee→WW full pipeline | Chiral vertices, triple gauge, massive propagators |

---

## WHAT EXISTS (v2, 329 tests)

### Six-layer pipeline

```
Layer 1: Model      → qed_model() / qcd_model() / ew_model()  → AbstractModel
Layer 2: Rules      → feynman_rules(model)                      → FeynmanRules
Layer 3: Channels   → tree_channels(model, rules, in, out)      → Vector{TreeChannel}
         Amplitude  → build_amplitude(ch, rules, model)          → DiracChains
Layer 4: Algebra    → trace → contract → expand → eval           → AlgSum (scalar)
Layer 5: Integrals  → PaVe{N}, evaluate(::PaVe; mu2)            → ComplexF64
Layer 6: Evaluate   → solve_tree(prob) → σ                      → Float64
```

### Pipeline coverage

| Process | Channels | Status |
|---------|----------|--------|
| e+e-→μ+μ- | s(γ) | Full pipeline via solve_tree |
| Bhabha | s(γ)+t(γ) | Pipeline + spin_sum_interference |
| Compton | s(e)+u(e) | Pipeline + _cross_line_trace |
| qq̄→gg | t(q)+u(q)+s(g) | Pipeline t/u, manual s (ggg vertex) |
| ee→WW | s(γ)+s(Z)+t(ν)+u(ν) | Channel enumeration only |

### Source files (35 files, ~3,600 LOC)

| File | LOC | What |
|------|-----|------|
| **Layer 4: Algebra** | | |
| `coeff.jl` | 142 | DimPoly coefficient algebra |
| `types.jl` | 79 | LorentzIndex, Momentum, MomentumSum |
| `colour_types.jl` | 126 | SUNT, SUNF, SUND, deltas |
| `pair.jl` | 76 | Parametric Pair{A,B} |
| `expr.jl` | 149 | AlgSum (Dict), AlgFactor, FactorKey |
| `sp_context.jl` | 71 | SPContext + ScopedValues |
| `contract.jl` | 136 | Lorentz contraction + Eps handling |
| `eps_contract.jl` | 85 | ε·ε = -det[pair(aᵢ,bⱼ)] |
| `expand_sp.jl` | 82 | Scalar product bilinear expansion |
| `dirac.jl` | 118 | DiracGamma{S}, Spinor{K}, DiracChain |
| `dirac_trace.jl` | ~160 | Trace → AlgSum (gamma5, projectors) |
| `dirac_expr.jl` | 104 | DiracExpr: matrix-valued expressions |
| `dirac_trick.jl` | 117 | D-dim γ^μ...γ_μ for n=0..5+ |
| `spin_sum.jl` | 138 | Fermion spin sums (completeness) |
| `interference.jl` | 98 | Cross-line traces, spin_sum_interference |
| `colour_trace.jl` | 72 | SU(N) trace → (real, imag) AlgSum |
| `colour_simplify.jl` | 148 | Delta contraction via dispatch |
| `polarization_sum.jl` | 58 | Feynman/axial/massive pol sums |
| **Layers 1-3** | | |
| `model.jl` | 99 | AbstractModel, QEDModel, Field{Species} |
| `qcd_model.jl` | 60 | QCDModel, qqg + ggg vertices |
| `ew_model.jl` | 49 | EWModel, 5 SM vertex types |
| `rules.jl` | 82 | FeynmanRules callable, vertex dispatch |
| `diagrams.jl` | 22 | ExternalLeg (mass field) |
| `channels.jl` | 103 | TreeChannel, tree_channels() |
| `amplitude.jl` | 141 | build_amplitude: boson + fermion exchange |
| **Layer 5: Integrals** | | |
| `pave.jl` | ~80 | PaVe{N} type, named constructors |
| `pave_eval.jl` | ~200 | evaluate: A₀/B₀/B₁ + C₁/C₂ PV reduction + _C0_quadgk |
| `c0_analytical.jl` | 82 | **C₀: COLLIER ccall + C0p0 analytical + quadgk fallback** |
| **Layer 6 + Reference** | | |
| `cross_section.jl` | ~100 | Mandelstam, solve_tree, σ |
| `schwinger.jl` | ~50 | REFERENCE: Schwinger correction |
| `vertex.jl` | 69 | REFERENCE: QED g-2 F₂(0)=α/(2π) |
| `running_alpha.jl` | ~100 | REFERENCE: running α(q²) |
| `ew_parameters.jl` | 33 | EW constants: M_W, M_Z, sin²θ_W |
| `ew_cross_section.jl` | 83 | REFERENCE: σ(ee→WW) Grozin formula |

### Test files (329 tests across 17 files)

| File | Tests | What |
|------|-------|------|
| `test_coeff.jl` | 29 | DimPoly arithmetic |
| `test_colour.jl` | 27 | SU(N) traces, δ contraction |
| `test_ee_mumu_x.jl` | 14 | e+e-→μ+μ- algebra (P&S 5.10) |
| `test_self_energy.jl` | 25 | DiracExpr, DiracTrick n=0,1,2 |
| `test_vertical.jl` | 34 | Full pipeline via solve_tree |
| `test_pave.jl` | 53 | PaVe types + A₀/B₀/B₁ |
| `test_schwinger.jl` | 15 | Schwinger correction + vacuum polarization |
| `test_compton.jl` | 4 | Compton |M|² from pipeline vs P&S 5.87 |
| `test_munit_batch1.jl` | 23 | MUnit: DiracTrace, Contract, PolarizationSum |
| `test_munit_batch2.jl` | 18 | MUnit: DiracTrick n=3,4 |
| `test_bhabha.jl` | 4 | Bhabha |M̄|² |
| `test_qqbar_gg.jl` | 2 | QCD qq̄→gg |M̄|² |
| `test_self_energy_1loop.jl` | 13 | 1-loop Σ(p) |
| `test_vertex_g2.jl` | 32 | C₀/C₁/C₂, F₂(0)=α/(2π) (SLOW: vertex_f2 uses own quadgk) |
| `test_running_alpha.jl` | 34 | Running α(q²), Δα, improved Born σ |
| `test_ee_ww.jl` | 36 | Tree-level e⁺e⁻→W⁺W⁻ reference formula |
| `test_pipeline.jl` | 17 | Bhabha/Compton/qq→gg/ee→WW pipeline |

---

## REFERENCE CODEBASES

All in `refs/` (gitignored):
- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — FF library source, Denner C0func.F backup.
- `refs/COLLIER/COLLIER-1.2.8/` — **COLLIER library. Build libcollier.so here.**
- `refs/papers/` — 17+ local paper copies incl. vanOldenborgh1990.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show master

# Build COLLIER (required per machine)
cd refs/COLLIER/COLLIER-1.2.8 && mkdir -p build && cd build
cmake .. -DCMAKE_Fortran_COMPILER=gfortran && make -j$(nproc)
cd /path/to/Feynfeld.jl && mkdir -p output

# Run specific test (fast)
julia --project=. test/v2/test_vertical.jl    # 5s, pipeline
julia --project=. test/v2/test_pipeline.jl     # 7s, all processes
julia --project=. test/v2/test_pave.jl         # 5s, integrals

# Fast smoke test (~40s)
for f in test/v2/test_coeff.jl test/v2/test_colour.jl test/v2/test_vertical.jl \
         test/v2/test_pipeline.jl test/v2/test_pave.jl; do
    julia --project=. "$f"
done

# Full suite (~10 min, vertex_g2 slow due to vertex_f2 reference quadgk)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
