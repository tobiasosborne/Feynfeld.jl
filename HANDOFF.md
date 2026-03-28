# HANDOFF — 2026-03-28 (End of Session 7, Spirals 5-7 complete)

## DO NOT DELETE THIS FILE. Read it completely before working.

## START HERE

1. Read `CLAUDE.md` — rules, anti-hallucination citation pattern, Julia idiom cheatsheet
2. Read `Feynfeld_PRD.md` — vision, spiral methodology, MUnit coverage target
3. Read `src/v2/DESIGN.md` — v2 type system, anti-patterns, cockroaches found
4. Read `JULIA_PATTERNS.md` — Julia idiom cheatsheet (same content as in CLAUDE.md §6)
5. Run `bd ready` to see available work (if beads errors, run `bd init --force --prefix feynfeld && bd backup restore`)
6. Run `for f in test/v2/test_*.jl; do julia --project=. "$f"; done` to verify 301 tests pass
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
| 7 | EW tree-level e+e-→W+W- | DONE (Session 7) |
| **8** | **MUnit mop-up** | **NEXT** |
| 7 | EW e+e-→W+W- | Planned |
| 8 | MUnit mop-up | Planned |
| 9+ | BSM / ULDM | Planned |

### Branch and code location

- **Branch:** `experimental/rebuild-v2`
- **v2 source:** `src/v2/` (28 files, ~3,400 LOC)
- **v2 tests:** `test/v2/` (16 files, 301 tests)
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
| `test_ee_ww.jl` | 33 | Tree-level e⁺e⁻→W⁺W⁻ at LEP2 energies |

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

### Spiral 7 completed: Tree-level e⁺e⁻ → W⁺W⁻

1. **EW parameters** (`ew_parameters.jl`): M_W, M_Z, sin²θ_W, Z-electron couplings
   (g_V, g_A, g_L, g_R). PDG 2024 on-shell values.

2. **Grozin analytical formula** (`ew_cross_section.jl`): Total cross-section from
   3 diagrams (s-channel γ/Z + t-channel ν_e) in massless electron limit. Formula
   from Grozin "Using REDUCE in HEP" Ch. 5.4, cross-checked via FeynCalc AnelEl-WW.m.
   Results at LEP2 energies: σ(200 GeV) ≈ 18 pb (tree-level; LEP measured ~16-17 pb
   including NLO corrections).

3. **Massive polarization sum**: `-g^μν + k^μk^ν/M²` for external W/Z bosons.

4. **33 tests**: EW parameters, massive pol sum, threshold behavior, LEP2 energies,
   high-energy gauge cancellation, α-scaling, M_W-dependence, QED comparison.

---

## WHAT TO DO NEXT

### Priority 1: Cleanup from Spiral 5-7 reviewer findings

These are small, self-contained fixes flagged by reviewers. Do them first.

1. **Fix massive pol sum citation** (`src/v2/polarization_sum.jl`, line ~40):
   Current: "Peskin & Schroeder, Eq. (5.75) generalized to massive case"
   Problem: P&S 5.75 is the massless axial-gauge sum, not the massive one.
   Fix: Cite `refs/FeynCalc/Tests/Feynman/PolarizationSum.test` IDs 2-3.
   Add verbatim equation string per Rule 2.

2. **Strengthen massive pol sum test** (`test/v2/test_ee_ww.jl`, lines 38-48):
   Current test only checks `isa AlgSum` and `length(terms) == 2`.
   Fix: Add test verifying actual coefficient values (metric tensor = -1,
   k^μk^ν term = 1/M²) by evaluating at a concrete kinematic point.

3. **Acquire Grozin book** (or demote citation):
   `ew_cross_section.jl` cites Grozin "Using REDUCE in HEP" Ch. 5.4 but the
   book is not in `refs/papers/`. Either fetch it (if available on arXiv/publisher)
   or change citation to secondary, with FeynCalc AnelEl-WW.m as primary.

### Priority 2: Spiral 8 — MUnit mop-up

Systematic translation of FeynCalc MUnit tests for functions already implemented.
Goal: push from 301 tests toward 500+.

**Protocol per MUnit test file** (from CLAUDE.md):
1. Read the `.test` file in `refs/FeynCalc/Tests/`.
2. Translate each `Test[]` to a Julia `@test`, preserving math exactly.
3. Document source file and test ID in a comment.
4. Cite the textbook equation that validates the test (Rule 2).
5. If MUnit test and textbook disagree, the textbook wins (Rule 1).

**Priority functions for translation:**

| MUnit file | Location | Est. tests | Notes |
|------------|----------|-----------|-------|
| `DiracTrace.test` | `Tests/Dirac/` | 50+ | Core function, many edge cases |
| `Contract.test` | `Tests/Lorentz/` | 40+ | Core Lorentz contraction |
| `PolarizationSum.test` | `Tests/Feynman/` | 15+ | Including massive vectors |
| `DiracTrick.test` | `Tests/Dirac/` | 30+ | General n≥5 case (already have n=0..4) |
| `ToPaVe.test` | `Tests/LoopIntegrals/` | 20+ | PaVe conversion tests |
| `PaVeReduce.test` | `Tests/LoopIntegrals/` | 20+ | C-function reductions |

**Approach:** Create `test/v2/test_munit_batch3.jl`, `batch4.jl`, etc. Each batch
should target ~20-30 tests from one or two MUnit files. Keep each test file < 200 LOC.

### Priority 3: Spiral 9+ — future work

| Spiral | Process | Key capability needed |
|--------|---------|----------------------|
| 9 | Full EW Feynman rules | `EWModel`, automated WWγ/WWZ vertices, EW diagram generation |
| 10 | BSM / ULDM | New model types, dark photon mixing, scalar DM |
| 11 | D₀ evaluation | 4-point scalar integral (box diagrams) |
| 12 | Full NLO automation | Tensor reduction + renormalization pipeline |

### Known open bugs

- **`feynfeld-83m`** — colour_trace n≥4 missing i² factor. The recursive
  decomposition T^aT^b = δ^{ab}/(2N) + (1/2)(d+if)T^c omits the imaginary
  unit i on the f term. When two f terms combine, i²=-1 is missing. Workaround:
  use analytical colour factors (Casimir identities) instead of recursive trace.
  See Session 6 HANDOFF for details.

- **B₀ imaginary part** — `_B0_quadgk` returns real part correctly but imaginary
  part from iε prescription is approximate (uses `complex(f, -1e-30)` rather
  than proper analytic continuation). Real part verified at <1e-8 accuracy.

- **vertex_f2 iε** — `vertex_f2(q2, m2, lambda2)` doesn't implement iε for
  q² > 4m² (above pair threshold). Documented in docstring. Use only for
  q² ≤ 0 (spacelike) or q² < 4m² (below threshold).

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
| `Denner1993_FortschrPhys41.pdf` | One-loop techniques, PV reduction, self-energies, EW |
| `PassarinoVeltman1979_NuclPhysB160.pdf` | PV decomposition |
| `tHooftVeltman1979_NuclPhysB153.pdf` | Scalar integrals (A₀, B₀, C₀, D₀) |
| `Shtabovenko2016_FeynCalc9_1601.01167.pdf` | FeynCalc 9.0 |
| `Shtabovenko2020_FeynCalc93_2001.04407.pdf` | FeynCalc 9.3 |
| `Shtabovenko2024_FeynCalc10_2312.14089.pdf` | FeynCalc 10 |
| `PDG2024_rev_standard_model.pdf` | EW review: α, sin²θ_W, M_W/Z, Eqs. (10.11)-(10.63) |
| `PDG2024_rev_quark_masses.pdf` | Quark masses review (MS-bar) |
| `PDG2024_rev_phys_constants.pdf` | Physical constants table |
| `PDG2024_sum_leptons.pdf` | Lepton summary: m_e, m_μ, m_τ |
| `PDG2024_sum_quarks.pdf` | Quark summary: m_u through m_t |
| `PDG2024_list_electron.pdf` | Full electron listings |
| `PDG2024_list_muon.pdf` | Full muon listings |
| `PDG2024_list_tau.pdf` | Full tau listings |
| `PDG2024_list_z_boson.pdf` | Full Z boson listings |

### Ground truth acquisition

See `CLAUDE.md` §Ground truth acquisition. Priority: arXiv → TIB VPN → playwright-cli.
The P&S textbook is at: `refs/papers/(Frontiers_in_Physics)...Peskin...(1995).djvu`
Note: .djvu format, unreadable by tools. Use FeynCalc examples as citation bridge.

---

## QUICK COMMANDS

```bash
# Branch
git branch  # should show experimental/rebuild-v2

# Run all v2 tests (301 tests across 16 files)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# Run specific test
julia --project=. test/v2/test_ee_ww.jl

# Fast smoke test (skip slow C₀/vertex tests, ~2 min total)
for f in test/v2/test_coeff.jl test/v2/test_colour.jl test/v2/test_vertical.jl \
         test/v2/test_schwinger.jl test/v2/test_running_alpha.jl test/v2/test_ee_ww.jl; do
    julia --project=. "$f"
done

# Beads
bd ready              # available work
bd stats              # project health
bd list --status=open # all open issues
bd create --title="..." --description="..." --type=task --priority=1

# Ground truth
ls refs/FeynCalc/Tests/                              # MUnit test directories
ls refs/FeynCalc/FeynCalc/Examples/                  # FeynCalc examples (QED/QCD/EW)
ls refs/papers/                                      # local paper copies (16 files)

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
├── ew_parameters.jl      # EW SM: M_W, M_Z, sin²θ_W, Z couplings (PDG 2024)
├── ew_cross_section.jl   # σ(e+e-→W+W-) Grozin formula, tree-level
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
├── test_running_alpha.jl      # Running α, Δα, improved Born σ (34 tests)
└── test_ee_ww.jl              # e⁺e⁻→W⁺W⁻ at LEP2 energies (33 tests)

Total: ~3,400 source LOC, ~2,000 test LOC, 301 tests, all files < 200 LOC.
```
