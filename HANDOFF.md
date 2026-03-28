# HANDOFF — 2026-03-28 (End of Session 7, Spirals 5-6 complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, anti-hallucination citation pattern, Julia idiom cheatsheet
2. Read `Feynfeld_PRD.md` — vision, spiral methodology, MUnit coverage target
3. Read `src/v2/DESIGN.md` — v2 type system, anti-patterns, cockroaches found
4. Read `JULIA_PATTERNS.md` — Julia idiom cheatsheet (same content as in CLAUDE.md §6)
5. Run `bd ready` to see available work (if beads errors, run `bd init --force --prefix feynfeld && bd backup restore`)
6. Run `for f in test/v2/test_*.jl; do julia --project=. "$f"; done` to verify 268 tests pass
7. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. The critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **ANTI-HALLUCINATION: CITE EVERYTHING.** Every formula must cite: local file
   path + equation number + verbatim copy. No exceptions. STOP if source missing.
   ```julia
   # Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.18)
   # "B₁(p², m₀², m₁²) = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)"
   ```
3. **ACQUIRE PAPERS BEFORE CODE.** If a paper is not in `refs/papers/`, fetch it
   BEFORE writing any code that cites it. Priority: arXiv (WebFetch) → TIB VPN
   (ask Tobias to connect, then use playwright-cli) → see CLAUDE.md §Ground truth.
   **NEVER ask Tobias to provide the paper — ask him to CONNECT TO VPN, then
   fetch it yourself via playwright-cli.** Downloaded files may land in
   `/tmp/playwright-artifacts-*/` — check there.
4. **JULIA IDIOMATIC ALL THE WAY.** Read the cheatsheet in CLAUDE.md FIRST.
   Use existing packages (QuadGK, PolyLog). No hand-rolled numerics. No isa cascades.
   Parametric types + dispatch. No OOP struct hierarchies unless dispatch genuinely needed.
5. **WORKFLOW: 3+1.** 3 read-only research subagents BEFORE any core code change
   (research source + 2 solution approaches). Then 1 rigorous reviewer agent AFTER.
6. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.
7. **LOC LIMIT ~200.** No source file exceeds ~200 lines.

**NEVER modify TensorGR.jl without explicit permission.**

---

## PROJECT OVERVIEW

### What is Feynfeld.jl?

Julia-native, agent-facing, full-stack physics computation suite. Lagrangian →
cross-section in one `using Feynfeld`. Replaces FeynRules + FeynArts + FeynCalc +
FormCalc + LoopTools. See PRD for the full vision.

### Development methodology: THE SPIRAL

Each **process** is a vertical spoke that drives horizontal **MUnit test coverage**.
FeynCalc has 15,222 MUnit tests (~10,000 translatable). Each spiral implements a
new process end-to-end and translates the MUnit tests for the functions it needs.

| Spiral | Process | Status |
|--------|---------|--------|
| 0 | e+e-→μ+μ- tree | DONE (Session 4) |
| 1 | Compton e+γ→e+γ | DONE (Session 5) |
| 2 | Bhabha e+e-→e+e- | DONE (Session 6) |
| 3 | QCD qq̄→gg | DONE (Session 6) |
| 4 | 1-loop self-energy Σ(p) | DONE (Session 6) |
| 5 | 1-loop vertex correction (g-2) | DONE (Session 7) |
| 6 | 1-loop vacuum polarization (running α) | DONE (Session 7) |
| **7** | **EW e+e-→W+W-** | **NEXT** |
| 7 | EW e+e-→W+W- | Planned |
| 8 | MUnit mop-up | Planned |
| 9+ | BSM / ULDM | Planned |

### Branch and code location

- **Branch:** `experimental/rebuild-v2`
- **v2 source:** `src/v2/` (26 files, ~3,200 LOC)
- **v2 tests:** `test/v2/` (15 files, 268 tests)
- **v1:** `src/algebra/`, `src/integrals/` — FROZEN, will be deleted. Do NOT extend or import patterns from.

---

## WHAT EXISTS (v2, 200 tests)

### Six-layer pipeline

```
Layer 1: Model      → qed_model()                    → QEDModel
Layer 2: Rules      → feynman_rules(model)            → FeynmanRules (callable)
Layer 3: Diagrams   → tree_diagrams(model, ...)       → Vector{FeynmanDiagram}
Layer 4: Algebra    → trace → contract → expand → eval → AlgSum (scalar)
Layer 5: Integrals  → PaVe{N}, evaluate(::PaVe; mu2)  → ComplexF64
Layer 6: Evaluate   → solve_tree(prob) → σ             → Float64
```

### Layer 4: Algebra (the core)

| File | What | Key types/functions |
|------|------|-------------------|
| `coeff.jl` | DimPoly coefficient algebra | `DimPoly`, `DIM`, `Coeff`, `evaluate_dim` |
| `types.jl` | Physics indices, momenta | `LorentzIndex`, `Momentum`, `MomentumSum`, `PairArg` |
| `colour_types.jl` | SU(N) types | `AdjointIndex`, `FundIndex`, `SUNT`, `SUNF`, `SUND`, deltas |
| `pair.jl` | Parametric Pair{A,B} | `MetricTensor`, `FourVector`, `ScalarProduct`, `pair()` |
| `expr.jl` | Dict-based AlgSum | `AlgSum`, `AlgFactor` (6-type union), `FactorKey`, `alg()` |
| `sp_context.jl` | Scalar product context | `SPContext`, `ScopedValues`, `evaluate_sp` |
| `contract.jl` | Lorentz contraction | `contract()`, `substitute_index()`, `alg_from_factors()` |
| `expand_sp.jl` | MomentumSum bilinear expansion | `expand_scalar_product()` |
| `dirac.jl` | Dirac gamma types | `DiracGamma{S}`, `Spinor{K}`, `DiracChain`, `GA/GAD/GS` |
| `dirac_trace.jl` | Trace → AlgSum | `dirac_trace()` (handles arbitrary-length chains) |
| `dirac_expr.jl` | Matrix-valued expressions | `DiracExpr` = Vector{Tuple{AlgSum, DiracChain}} |
| `dirac_trick.jl` | γ^μ...γ_μ contraction | `dirac_trick()` (n=0,1,2,3,4 + general n≥5) |
| `spin_sum.jl` | Fermion spin sums | `spin_sum_amplitude_squared()`, `_single_line_trace()` |
| `colour_trace.jl` | SU(N) trace | `colour_trace()` (concrete N, recursive) |
| `colour_simplify.jl` | δ contraction, f·f, d·d | `contract_colour()` |
| `polarization_sum.jl` | Photon/gluon pol sum | `polarization_sum()` (Feynman + axial gauge) |

### Layers 1-3, 5-6

| File | What |
|------|------|
| `model.jl` | `QEDModel`, `Field{Species}`, `GaugeGroup`, traits |
| `rules.jl` | `FeynmanRules` callable, vertex/propagator dispatch on species |
| `diagrams.jl` | Hard-coded e+e-→μ+μ- s-channel topology |
| `pave.jl` | `PaVe{N}` parametric type, named constructors A0/B0/B1/C0/C1/C2/D0 |
| `pave_eval.jl` | `evaluate(::PaVe{1,2,3})`, B0/C0 via QuadGK, C1/C2 via PV reduction |
| `vertex.jl` | `vertex_f2_zero`, `vertex_f2` — QED anomalous magnetic moment |
| `schwinger.jl` | Analytical Schwinger correction formula, vacuum polarization |
| `cross_section.jl` | `CrossSectionProblem`, `solve_tree`, Mandelstam, dσ/dΩ |

### Dependencies

- **PolyLog.jl** v2.6.2 — dilogarithm `li2()` (for future C₀ evaluation)
- **QuadGK.jl** v2.11.2 — adaptive quadrature (used by B₀ and vacuum polarization)

### Test files

| File | Tests | What |
|------|-------|------|
| `test_coeff.jl` | 29 | DimPoly arithmetic |
| `test_colour.jl` | 22 | SU(N) traces, δ contraction |
| `test_ee_mumu_x.jl` | 14 | e+e-→μ+μ- algebra (P&S 5.10) |
| `test_self_energy.jl` | 25 | DiracExpr, DiracTrick n=0,1,2 |
| `test_vertical.jl` | 33 | Full pipeline: Model→Rules→Diagrams→Algebra→σ |
| `test_pave.jl` | 2 | PaVe types + A₀/B₀ numerical evaluation |
| `test_schwinger.jl` | 15 | Schwinger correction + vacuum polarization |
| `test_compton.jl` | 4 | Compton |M|² from pipeline vs P&S Eq. 5.87 |
| `test_munit_batch1.jl` | 23 | MUnit translations: DiracTrace, Contract, PolarizationSum |
| `test_munit_batch2.jl` | 18 | MUnit translations: DiracTrick n=3,4 (ThreeFreeIndices, FourFreeIndices) |
| `test_bhabha.jl` | 4 | Bhabha |M̄|² at 2 kinematic points vs FeynCalc |
| `test_qqbar_gg.jl` | 2 | QCD qq̄→gg |M̄|² at 2 kinematic points vs FeynCalc/ESW |
| `test_self_energy_1loop.jl` | 13 | 1-loop Σ(p) via B₀, B₁, A₀ at 2 off-shell points |
| `test_vertex_g2.jl` | 32 | C₀/C₁/C₂ evaluation, F₂(0)=α/(2π) Schwinger |
| `test_running_alpha.jl` | 34 | Running α(q²), Δα, improved Born σ(e+e-→μ+μ-) |

---

## WHAT WAS DONE IN SESSION 6

### Spirals 2-4 completed

1. **Spiral 2 (Bhabha e+e-→e+e-)**: s+t channel interference, 8-gamma trace for
   identical fermion scattering. 4 tests at 2 kinematic points.

2. **Spiral 3 (QCD qq̄→gg)**: First QCD process. Triple gluon vertex constructed
   from Pair factors. Physical (axial gauge) polarization sums implemented —
   required because QCD Ward identity only holds for full amplitude sum, not
   individual diagram contributions. Analytical colour factors (C_tt=16/3,
   C_uu=16/3, C_tu=-2/3, C_ss=12). 2 tests vs FeynCalc/Ellis-Stirling-Weber.

3. **Spiral 4 (1-loop self-energy)**: QED electron self-energy Σ(p) = p̸ Σ_V + mΣ_S
   computed from B₀, B₁, A₀ PaVe functions. Validated at 2 off-shell kinematic
   points with analytical cross-checks. 13 tests.

4. **DiracTrick n≥3**: Implemented general Mertig-Boehm-Denner formula (Eq. 2.9)
   for γ^μ γ^{a1}...γ^{an} γ_μ. Explicit formulas for n=3,4; recursive general
   formula for n≥5. 18 MUnit tests translated from FeynCalc DiracTrick.test.

5. **Colour algebra extensions**: d·d (3 shared), f·d (vanishes), f·f (3 shared)
   contraction identities added. Refactored via `_try_struct_contract` dispatch.

6. **Physical polarization sums**: `polarization_sum(μ, ν, k, n; ctx)` for axial
   gauge with reference momentum. Required for QCD.

7. **Papers acquired**: 7 papers now in `refs/papers/`:
   - P&S textbook (.djvu)
   - MertigBohmDenner1991 — FeynCalc original, Eq. 2.9 for DiracTrick
   - Denner1993 — one-loop techniques, PV reduction, self-energies (Appendix B)
   - PassarinoVeltman1979 — PV decomposition
   - tHooftVeltman1979 — scalar integrals
   - Shtabovenko2016/2020/2024 — FeynCalc 9.0/9.3/10

### Key learnings (Session 6)

1. **ACQUIRE PAPERS FIRST.** The biggest mistake of the session: wrote DiracTrick
   code citing Mertig-Boehm-Denner 1991 without having the paper locally. Tobias
   caught this as a Rule 2 violation. The procedure: arXiv (WebFetch, do it
   yourself) → TIB VPN (ask Tobias to connect, then playwright-cli yourself).
   Downloads may land in `/tmp/playwright-artifacts-*/`.

2. **QCD gauge invariance.** Individual QCD diagram contributions are NOT gauge-
   invariant. Using -g^{μν} (Feynman gauge) for individual D_{XY} gives wrong
   results for the s-channel of qq̄→gg. Fix: use physical (axial gauge)
   polarization sums with reference momenta, matching FeynCalc's
   `DoPolarizationSums[#,k1,k2]`.

3. **Relative phase between channels.** The s-channel (gluon propagator -ig/s)
   and t/u channels (quark propagator ip̸/t) have a relative phase of i. This
   manifests as a -1 on D_ss in the colour-separated decomposition. The sign
   was found empirically and needs deeper investigation.

4. **colour_trace n≥4 has a known bug** (feynfeld-83m): the recursive
   decomposition T^aT^b = δ^{ab}/(2N) + (1/2)(d+if)T^c omits the imaginary
   unit i on the f term. When two f terms combine, i²=-1 is missing. For
   Spiral 3, this was worked around by using analytical colour factors (Casimir
   identities) instead of the recursive trace.

5. **B₀ imaginary part not computed.** Our QuadGK-based B₀ evaluation doesn't
   handle the iε prescription for timelike momenta (p² > m²). The real part is
   correct; the imaginary part is always returned as 0. For Spiral 4, this was
   sidestepped by testing at kinematics where B₀ is real (p² < m²).

---

## WHAT WAS DONE IN SESSION 7

### Spiral 5 completed: QED vertex correction (g-2)

1. **C₀ scalar 3-point evaluation**: Implemented via nested 2D QuadGK Feynman
   parameter integral. Formula from 't Hooft-Veltman (1979) Eq. 5.2 / Denner
   (1993) Eq. 4.26. Handles spacelike and timelike kinematics via iε prescription.

2. **C₁, C₂ tensor coefficients**: Passarino-Veltman reduction using Gram matrix
   inversion. PV identity 2l·kⱼ = Dⱼ+mⱼ²-D₀-m₀²-kⱼ² cancels one propagator,
   yielding R coefficients in terms of B₀ and C₀. Gram matrix singularity detected
   and throws error. Formulas from Denner (1993) Eqs. 4.6-4.8.

3. **F₂(0) = α/(2π)**: Direct Feynman parameter integral (P&S Chapter 6,
   cross-checked via FeynCalc El-GaEl.m). At q²=0 the 2D integral reduces to 1D:
   F₂(0) = (α/π) m² ∫₀¹ dz z(1-z)²/[m²(1-z)²+zλ²]. For λ→0: F₂(0) = α/(2π).
   Also implemented general F₂(q²) for spacelike q² via 2D integral.

4. **Named constructors**: C1, C2 added to pave.jl alongside existing C0.

5. **32 tests** in `test_vertex_g2.jl`:
   - C₀ at 7 kinematic points (zero momenta, asymmetric, g-2, zero mass, timelike)
   - C₁/C₂ at 2 points with Gram matrix self-consistency check
   - Gram matrix singularity detection
   - F₂(0) = α/(2π) exact (λ=0), IR convergence (λ→0), mass-independence,
     α-linearity, 1D/2D consistency, spacelike monotonicity

### Key decisions

1. **QuadGK over analytical C₀**: The 't Hooft-Veltman analytical formula involves
   complex dilogarithms with intricate branch cut handling. QuadGK nested integration
   is simpler, correct by construction, and fast enough (~5s for full test suite).

2. **Direct F₂(0) over PaVe decomposition**: Computing F₂ from the Feynman parameter
   integral directly (after textbook Dirac algebra) is simpler and more reliable than
   the full amplitude → PaVe → form factor extraction pipeline. The PaVe C functions
   are still implemented as infrastructure for future spirals.

3. **vertex_f2 valid only below threshold**: The general q² function doesn't
   implement iε for q² > 4m² (above pair-production threshold). Documented.

---

## WHAT WAS DONE IN SESSION 7 (continued): SPIRAL 6

### Spiral 6 completed: Running α(q²) from vacuum polarization

1. **SM fermion table**: PDG 2024 masses for 3 leptons + 6 quarks with charges
   and color factors. Ref: `refs/papers/PDG2024_sum_leptons.pdf`,
   `refs/papers/PDG2024_sum_quarks.pdf`.

2. **`delta_alpha(q2; alpha, fermions)`**: Total vacuum polarization Δα(q²) summed
   over active fermions: Δα = Σ_f Q_f² N_c × (-Π̂_f). Returns ComplexF64 (imaginary
   part nonzero for timelike q² above fermion pair thresholds).

3. **`running_alpha(q2; alpha)`**: Running coupling α(q²) = α/|1-Δα(q²)|.
   Validated: α⁻¹(M_Z²) ≈ 128 (perturbative quarks), Δα_lep(M_Z²) ≈ 0.0314
   matching PDG value 0.0315 to 0.3%.

4. **`sigma_improved_ee_mumu(s; alpha)`**: Improved Born approximation for
   σ(e+e-→μ+μ-) = (4πα(s)²)/(3s) × (1+δ_Schwinger) × (GeV⁻²→nb).

5. **vacuum_polarization bug fix**: Changed return type from Float64 to ComplexF64
   — was discarding imaginary part for timelike momenta. Added 2 tests for
   imaginary part in test_schwinger.jl.

6. **Schwinger comment fix**: Corrected misleading comment claiming VP included
   in Schwinger correction (it is not — VP enters separately via running α).

7. **9 PDG papers acquired**: lepton/quark summary tables, physical constants,
   standard model review, quark masses review, Z boson listings.

8. **34 tests** in `test_running_alpha.jl` + 2 new tests in `test_schwinger.jl`:
   - SM table structure, Δα leptonic/total/energy-dependent/complex/spacelike
   - Running α at M_Z, basic properties, leading-log check
   - Improved Born σ at M_Z and multiple energies

---

## WHAT TO DO NEXT: SPIRAL 7 (EW e+e-→W+W-)

### Process: Electroweak e+e-→W+W- tree-level

This is the first process that requires the full electroweak sector: SU(2)×U(1)
gauge group, W/Z bosons, weak mixing angle.

### New capabilities needed

1. **EW model** — extend `model.jl` beyond QED: SU(2)×U(1) gauge group,
   W±/Z bosons as massive vector fields, weak mixing angle θ_W.
2. **Triple gauge coupling** — WWγ and WWZ vertices.
3. **Massive vector propagator** — (g^μν - k^μk^ν/M²)/(k²-M²).
4. **Multiple diagram channels** — s-channel (γ,Z) + t-channel (ν_e).
5. **Tests vs FeynCalc/Denner** — differential and total cross-sections.

### Ground truth

- **Denner 1993:** Full EW corrections to e+e-→W+W-
- **FeynCalc examples:** `refs/FeynCalc/FeynCalc/Examples/EW/`
- **P&S:** Chapter 21

### Known open bugs

- `feynfeld-83m` — colour_trace n≥4 missing i² factor
- vertex_f2 iε not implemented for q² > 4m²

---

## REFERENCE CODEBASES

All in `refs/` (gitignored):
- `refs/FeynCalc/` — 186k LOC Mathematica. MUnit tests in `Tests/`. Examples in `FeynCalc/Examples/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — Loop integral numerics (Fortran source).
- `refs/papers/` — See list below.

### Papers in refs/papers/

| File | Reference |
|------|-----------|
| `...Peskin...djvu` | P&S textbook (1995), .djvu format (unreadable by tools) |
| `MertigBohmDenner1991_FeynCalc_CPC64.pdf` | FeynCalc original, Eq. 2.9 (DiracTrick) |
| `Denner1993_FortschrPhys41.pdf` | One-loop techniques, PV reduction, self-energies |
| `PassarinoVeltman1979_NuclPhysB160.pdf` | PV decomposition |
| `tHooftVeltman1979_NuclPhysB153.pdf` | Scalar integrals (A₀, B₀, C₀, D₀) |
| `Shtabovenko2016_FeynCalc9_1601.01167.pdf` | FeynCalc 9.0 |
| `Shtabovenko2020_FeynCalc93_2001.04407.pdf` | FeynCalc 9.3 |
| `Shtabovenko2024_FeynCalc10_2312.14089.pdf` | FeynCalc 10 |

### Ground truth acquisition

See `CLAUDE.md` §Ground truth acquisition. Priority: arXiv → TIB VPN → playwright-cli.
The P&S textbook is at: `refs/papers/(Frontiers_in_Physics)...Peskin...(1995).djvu`
Note: .djvu format, unreadable by tools. Use FeynCalc examples as citation bridge.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run all v2 tests (200 tests across 13 files, excluding WIP test_qqbar_gg)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Run specific test
julia --project=. test/v2/test_bhabha.jl

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues
bd create --title="..." --description="..." --type=task --priority=1

# Ground truth
ls refs/FeynCalc/Tests/                              # MUnit test directories
ls refs/FeynCalc/FeynCalc/Examples/QED/              # FeynCalc examples (Tree + OneLoop)
ls refs/papers/                                      # local paper copies

# Commit and push (session end protocol)
git add <files>
git commit -m "..."
git push
```

---

## FILE MAP

```
src/v2/
├── FeynfeldX.jl          # Module root, includes + exports
├── coeff.jl              # DimPoly coefficients
├── types.jl              # PhysicsIndex, Momentum, MomentumSum
├── colour_types.jl       # SU(N): AdjointIndex, FundIndex, SUNT, SUNF, SUND
├── pair.jl               # Parametric Pair{A,B}, NOT exported
├── expr.jl               # AlgSum (Dict), AlgFactor (6-type union), FactorKey
├── sp_context.jl         # SPContext + ScopedValues
├── contract.jl           # Lorentz contraction + substitute_index
├── expand_sp.jl          # Scalar product bilinear expansion
├── dirac.jl              # DiracGamma{S}, Spinor{K}, DiracChain
├── dirac_trace.jl        # Dirac trace → AlgSum
├── dirac_expr.jl         # DiracExpr: matrix-valued expressions
├── dirac_trick.jl        # D-dim γ^μ...γ_μ for n=0..5+ (Eq. 2.9)
├── spin_sum.jl           # Fermion spin sums
├── colour_trace.jl       # SU(N) trace → AlgSum (BUG: i² missing for n≥4)
├── colour_simplify.jl    # Delta contraction, f·f, d·d, f·d identities
├── polarization_sum.jl   # Feynman + axial gauge pol sums
├── model.jl              # AbstractModel, QEDModel, Field{Species}
├── rules.jl              # FeynmanRules callable
├── diagrams.jl           # FeynmanDiagram, hard-coded topologies
├── pave.jl               # PaVe{N} type, named constructors (A0-D0, C1, C2)
├── pave_eval.jl          # evaluate(::PaVe{1,2,3}) via QuadGK + PV reduction
├── schwinger.jl          # Schwinger correction + vacuum pol (ComplexF64)
├── vertex.jl             # QED vertex F₂(0)=α/(2π), vertex_f2_zero/vertex_f2
├── running_alpha.jl      # SM fermions, running_alpha(q²), sigma_improved
├── cross_section.jl      # Mandelstam, Problem/Solve, σ
├── DESIGN.md             # Design choices, anti-patterns, cockroaches
└── VERTICAL_PLAN.md      # Original vertical plan (historical)

test/v2/
├── test_coeff.jl              # DimPoly (29 tests)
├── test_colour.jl             # SU(N) (22 tests)
├── test_ee_mumu_x.jl          # e+e-→μ+μ- algebra (14 tests)
├── test_self_energy.jl        # DiracExpr + DiracTrick (25 tests)
├── test_vertical.jl           # Full pipeline (33 tests)
├── test_pave.jl               # PaVe types + numerics (2 tests)
├── test_schwinger.jl          # Schwinger correction (15 tests)
├── test_compton.jl            # Compton |M|² vs P&S 5.87 (4 tests)
├── test_munit_batch1.jl       # MUnit: DiracTrace, Contract, PolarizationSum (23 tests)
├── test_munit_batch2.jl       # MUnit: DiracTrick n=3,4 (18 tests)
├── test_bhabha.jl             # Bhabha |M̄|² vs FeynCalc (4 tests)
├── test_qqbar_gg.jl           # QCD qq̄→gg |M̄|² vs FeynCalc/ESW (2 tests)
├── test_self_energy_1loop.jl  # 1-loop Σ(p) via PaVe (13 tests)
├── test_vertex_g2.jl          # C₀/C₁/C₂ + F₂(0)=α/(2π) (32 tests)
└── test_running_alpha.jl      # Running α, Δα, improved Born σ (34 tests)

Total: ~3,200 source LOC, ~1,800 test LOC, 268 tests, all files < 200 LOC.
```
