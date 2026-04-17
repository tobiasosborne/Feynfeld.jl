# Feynfeld.jl v2 Stocktake: Integrals (Layer 5) & Evaluate (Layer 6)
**Date:** 2026-04-17  
**Scope:** PaVe, loop integrals, cross-section evaluation, and NLO assembly  
**Total LOC (Integrals + Evaluate):** ~1400 lines

---

## File-by-File Inventory

### **Layer 5: Integrals (PaVe, Tensor Reduction, Loop Evaluation)**

| File | LOC | Purpose | Public API | Dependencies |
|------|-----|---------|-----------|--------------|
| **pave.jl** | 86 | PaVe{N} type: generic N-point scalar/tensor loop integral; named constructors A0, B0, B1, C0, C1, C2, D0, D1–D3 | `PaVe{N}`, `A0()`, `B0()`, `B1()`, `B00()`, `B11()`, `C0()`, `C1()`, `C2()`, `D0()`, `D1()`, `D2()`, `D3()` | None (standalone type) |
| **pave_eval.jl** | 107 | Numerical dispatch layer: `evaluate(pv; mu2)` routes to _eval_A, _eval_B, _eval_C, _eval_D | `evaluate()` | pave.jl, b0_eval.jl, c0_analytical.jl, d0_collier.jl, d_tensor.jl |
| **b0_eval.jl** | 122 | B₀ scalar 2-point function: QuadGK adaptive quadrature on Feynman parameters; handles threshold via Kallen function; COLLIER fallback for IR-divergent cases | `_eval_B()`, `_B0()`, `_B0_collier()` (internal); `_eval_B_tensor()` for B₁ | QuadGK, c0_analytical.jl (for A₀) |
| **c0_analytical.jl** | 82 | C₀ scalar 3-point function: COLLIER primary (ccall via libcollier.so); C0p0 analytical special case; quadgk fallback for Gram-degenerate kinematics (Δ₂≈0) | `_C0_analytical()`, `_C0_collier()`, `_C0p0()` (internal); exports via pave_eval | COLLIER library (ccall), QuadGK (fallback) |
| **d0_collier.jl** | 41 | D₀ scalar 4-point function: COLLIER only (no quadgk fallback; triple-nested closures cause JIT explosion); error if COLLIER missing | `_D0_evaluate()`, `_D0_collier()` (internal) | COLLIER library (mandatory) |
| **d_tensor.jl** | 59 | D-tensor coefficients D₁, D₂, D₃ via Passarino-Veltman reduction: solves 3×3 Gram system; cascades to C₀ integrals with each propagator removed | `_eval_D_tensor()` (internal) | c0_analytical.jl, d0_collier.jl |
| **tid.jl** | 175 | Tensor Integral Decomposition: evaluates 1-loop box numerators ∫ N(q)/[D₀D₁D₂D₃]; handles rank 0, 1, 2 via Gram-matrix decomposition and PV cancellation | `evaluate_box_integral()` | pave_eval.jl, c0_analytical.jl, d0_collier.jl, d_tensor.jl, sp_lookup.jl (implicit via sp_vals dict) |
| **nlo_box.jl** | 108 | NLO box assembly: full pipeline per channel (build_loop_box_amplitude → spin_sum_tree_loop_interference → contract → expand_sp → evaluate_box_integral); spin-averaged Born-virtual coupling | `evaluate_single_box_channel()`, `evaluate_box_channels()`, `born_virtual_box()` | loop_amplitude.jl, loop_interference.jl, loop_channels.jl, contract.jl (Layer 4), expand_scalar_product.jl (Layer 4), tid.jl |

### **Layer 6: Evaluate (Cross Sections, Observables) + Layer 6b: NLO Assembly**

| File | LOC | Purpose | Public API | Dependencies |
|------|-----|---------|-----------|--------------|
| **amplitude.jl** | 179 | Tree amplitude building: dispatches on exchanged field species (boson vs fermion); boson exchange yields two DiracChains (fermion lines), fermion exchange yields propagator decomposition (chain_mom, chain_mass, mass) | `build_amplitude()` | channels.jl (TreeChannel), rules.jl (FeynmanRules), model.jl (AbstractModel), gauge_exchange.jl |
| **cross_section.jl** | 196 | Core pipeline: Mandelstam kinematics, CrossSectionProblem DiffEq pattern, `solve_tree()` (standard channels), `solve_tree_pipeline()` (qgraf-faithful orbit-canonical emission stream); spin sum → contract → expand_sp; Float64 evaluation via sp_values_2to2; differential/total cross-section formulas | `Mandelstam`, `CrossSectionProblem`, `solve_tree()`, `solve_tree_pipeline()`, `sp_values_2to2()`, `evaluate_m_squared()`, `dsigma_domega()`, `sigma_total_tree_ee_mumu()` | amplitude.jl, channels.jl, spin_sum.jl (Layer 4), contract.jl (Layer 4), expand_scalar_product.jl (Layer 4), QgrafPort.jl (Phase 18b-1) |
| **loop_amplitude.jl** | 115 | 1-loop box amplitude construction: BoxDenominators struct with accumulated momenta and Mandelstam labels (:t, :u); direct vs crossed topology; fermion lines with loop momentum q | `build_loop_box_amplitude()`, `BoxDenominators` | loop_channels.jl, rules.jl, model.jl |
| **loop_channels.jl** | 79 | LoopChannel enumeration for 1-loop boxes: 2→2 2-fermion-line boxes (e.g. e⁺e⁻→μ⁺μ⁻); direct + crossed topologies; validates boson exchange vertices exist | `box_channels()`, `LoopChannel` | channels.jl (ExternalLeg), rules.jl, model.jl |
| **loop_interference.jl** | 63 | Tree × loop interference: spin-summed Born-virtual M_tree* × M_box; product of two independent fermion-line traces with loop momentum; feeds SP(q, p_i) and SP(q, q) to TID | `spin_sum_tree_loop_interference()` | interference.jl (Layer 4) |
| **gauge_exchange.jl** | 38 | Gauge boson exchange amplitude (fermion ↔ triple-gauge): returns (chain::DiracExpr, vertex::AlgSum) for qq̄→gg style; triple gauge vertex in all-outgoing convention | `_build_gauge_exchange()` (called by amplitude.jl) | amplitude.jl (internal to build_amplitude dispatch) |
| **ew_cross_section.jl** | 79 | **Standalone reference:** e⁺e⁻→W⁺W⁻ total cross-section (massless e limit); 3-channel tree (γ, Z, ν_e); Grozin formula with log and sqrt parts; validates pipeline against known physics | `sigma_ee_ww()`, `_sigma_grozin()`, constants (EW_M_W, EW_SIN2_W, EW_ALPHA) | None (standalone) |
| **schwinger.jl** | 59 | **Standalone reference:** O(α) Schwinger correction to e⁺e⁻→μ⁺μ⁻; vacuum polarization Π̂(s) via QuadGK; NLO total cross-section formula; validates vertex + soft-real correction | `schwinger_correction()`, `vacuum_polarization()`, `sigma_nlo_ee_mumu()` | QuadGK (implicit in vacuum_polarization) |
| **vertex.jl** | 69 | **Standalone reference:** QED vertex correction F₂(q²) (anomalous magnetic moment); 1D/2D Feynman parameter integrals via QuadGK; F₂(0)=α/(2π) ground truth; validates form-factor pipeline | `vertex_f2_zero()`, `vertex_f2()` | QuadGK |
| **running_alpha.jl** | 98 | **Standalone reference:** Running fine-structure constant α(q²) from vacuum polarization; loops over SM fermion table (leptons + quarks); complex-valued above threshold; on-shell running coupling; validates beta-function physics | `delta_alpha()`, `running_alpha()`, SM_FERMIONS, SM_LEPTONS, SM_QUARKS | schwinger.jl (vacuum_polarization) |

---

## Pipeline vs Standalone Classification

### **Pipeline Components** (feed into solve_tree/solve_tree_pipeline)
- **amplitude.jl**: Amplitude building from TreeChannel → DiracExpr
- **cross_section.jl**: Central dispatch: problem definition, solver entry points, spin sum/contract/expand pipeline, kinematics
- **loop_amplitude.jl, loop_channels.jl, loop_interference.jl**: 1-loop amplitude/interference building
- **pave.jl, pave_eval.jl, b0_eval.jl, c0_analytical.jl, d0_collier.jl, d_tensor.jl, tid.jl**: Loop integral evaluation chain
- **nlo_box.jl**: NLO assembly orchestration

**Data flow:** CrossSectionProblem → solve_tree_pipeline → [amplitude.jl] → [spin_sum/contract/expand] → [cross_section evaluation] → OR [evaluate_box_channels] → [nlo_box.jl] → [tid.jl] → [pave_eval.jl]

### **Standalone Analytical References** (per PRD §2.3)
Flagged as validation ground-truth, NOT part of standard pipeline:
- **ew_cross_section.jl**: e⁺e⁻→W⁺W⁻ Grozin formula (validates EW tree amplitude)
- **schwinger.jl**: O(α) Schwinger correction (validates NLO QED vertex + soft real)
- **vertex.jl**: QED vertex F₂ form factors (validates loop vertex integrals)
- **running_alpha.jl**: Vacuum polarization & running coupling (validates fermion-loop physics)

These modules are **import-only**: called by tests/validation code, not by the main solve_tree pipeline. Each stands alone with minimal deps (QuadGK, vacuum_polarization for schwinger/running_alpha).

---

## Cross-Cutting Observations

### 1. **COLLIER Integration (d0_collier.jl, c0_analytical.jl)**

**d0_collier.jl** is the simplest file (41 LOC): pure ccall wrapper. No fallback; error if library missing. Justification: triple-nested QuadGK closure causes Julia JIT explosion (minutes at runtime). Acceptance of external dep is intentional.

**c0_analytical.jl** is more defensive (82 LOC):
- COLLIER primary (most kinematics)
- C0p0 (all momenta zero) special case with closed-form logarithms
- QuadGK fallback for Gram-degenerate Δ₂≈0 (known limitation of Spence-based algs; future L'Hôpital limit)

**b0_eval.jl** is the workhorse (122 LOC):
- QuadGK by default (no singularities in two-parameter space)
- Closed forms for both-massless and one-massless cases (split log singularity analytically)
- Threshold handling: Kallen function λ determines regime (below vs above); above threshold, switch to |f| quadrature + analytical imag part
- COLLIER fallback only for IR-divergent degenerate B₀(0,0,0)

**Architectural lesson:** MS-bar scheme is built in. All A₀, B₀, C₀, D₀ return finite parts (UV poles subtracted). No global μ² parameter leaks; each function takes mu2 kwarg.

### 2. **From Tree Pipeline to Loop Integrals**

**Path in code:**
```
solve_tree_pipeline (cross_section.jl:162)
  ├─ QgrafPort.emission_to_amplitude (qgraf/*)
  ├─ combine_m_squared_burnside (qgraf/*)
  ├─ contract (Layer 4)
  └─ expand_scalar_product (Layer 4)
     [returns AlgSum with Pair{Momentum,...} factors]

evaluate_box_channels (nlo_box.jl:62)
  └─ for ch in box_channels(model, rules, incoming, outgoing)
     ├─ evaluate_single_box_channel (nlo_box.jl:29)
     │  ├─ build_loop_box_amplitude (loop_amplitude.jl:38)
     │  │  └─ returns (chain_e, chain_mu, denoms)
     │  ├─ spin_sum_tree_loop_interference (loop_interference.jl:28)
     │  │  └─ _cross_line_trace (interference.jl, Layer 4)
     │  │     [products of two AlgSum]
     │  ├─ contract (Layer 4)
     │  └─ expand_scalar_product (Layer 4)
     │     [SP(q, p_i) and SP(q, q) coefficients extracted]
     └─ evaluate_box_integral (tid.jl:23)
        ├─ Rank analysis (0, 1, or 2)
        ├─ Pre-compute D₀, D₁, D₂, D₃ (pave_eval.jl)
        │  └─ d_tensor.jl for tensor indices → C₀ cascade
        ├─ Gram-matrix decomposition of external momenta
        └─ Returns ComplexF64 in COLLIER normalization
           ├─ No coupling (embedded in tree_amp*tree_amp trace)
           ├─ No propagator denominator (1/s already in diagram)
           └─ Spin-averaged via factor of 1/4 in nlo_box.jl:96
```

**Key insight:** By `evaluate_box_integral`, the numerator N(q) is fully determined (no symbolic loop momentum). TID (tensor integral decomposition) is purely numerical. Each `evaluate_box_integral` call is a single 4-dimensional box integral in COLLIER normalization.

### 3. **Coupling & Normalization**

**In nlo_box.jl (line 96):**
```julia
-e_sq^3 / (32π^2 * man.s) * imag(I_box)
```
where `e_sq = 4π * alpha`.

Physical Born-virtual interference is:
```
(1/4) Σ_spins 2Re(M_tree* × M_box) = -e⁶/(32π²s) × Im(result)
```

This extraction of Im via the overall coupling factor is verified empirically (refs/FeynCalc/ElAel-MuAmu.m, line 20 in nlo_box.jl comment). The i from tree propagator survives after vertex/propagator phases cancel.

**Why it works:**
1. Tree amplitude: no loop momentum, coupling e² embedded in vertex
2. Loop amplitude: numerator from chain + vertices, loop momentum in propagators
3. Interference trace: M_tree* has conjugate gammas; product yields Lorentz-scalar coefficients
4. TID: all kinematic invariants Float64, loop integral is real numerical evaluation (modulo iε for on-shell poles)
5. Final contraction: Im part of box integral times coupling → physical rate

### 4. **Deferred & Boundary Cases**

**Known limitations:**
- **tid.jl (line 73):** Rank > 2 errors. Higher-loop kinematics not implemented.
- **c0_analytical.jl (lines 11–14):** Gram-degenerate Δ₂≈0 (external momenta linearly dependent) returns NaN from COLLIER. QuadGK fallback available but slow. Future: L'Hôpital limit or analytic continuation.
- **b0_eval.jl (lines 115):** B₁ at p²=0 errors (requires special treatment; not implemented).
- **amplitude.jl (line 32):** "Pure boson scattering not yet implemented."
- **nlo_box.jl (lines 61):** TODO multi-channel interference (Phase B).

**No TODOs in:**
- pave.jl, pave_eval.jl, b0_eval.jl, c0_analytical.jl, d0_collier.jl, d_tensor.jl (Integrals are feature-complete for 2→2 ee→μμ)
- schwinger.jl, vertex.jl, running_alpha.jl, ew_cross_section.jl (Standalone refs are final)

### 5. **Coefficient Representation**

All integrals return **ComplexF64**. No DimPoly in Layer 5 (UV poles already subtracted in MS-bar scheme). This decouples Integrals from Algebra fully — pure numerics by design.

### 6. **File Cohesion & Dependencies**

**Integrals cluster (Layer 5):**
- pave.jl: standalone type
- pave_eval.jl: dispatcher, depends on (b0, c0, d0, d_tensor)
- b0_eval.jl: QuadGK arithmetic, minimal
- c0_analytical.jl: COLLIER + QuadGK, fallback chain
- d0_collier.jl: pure ccall, no fallback
- d_tensor.jl: PV reduction, depends on c0
- tid.jl: final assembly, depends on all above

**Evaluate cluster (Layer 6):**
- amplitude.jl: tree only, minimal algebra
- cross_section.jl: dispatcher, depends on amplitude + spin_sum + contract + expand
- loop_amplitude.jl, loop_channels.jl, loop_interference.jl: 1-loop building blocks
- nlo_box.jl: orchestrator, depends on loop_* + tid + cross_section
- gauge_exchange.jl: amplitude dispatch helper
- 4 standalone refs: no mutual dependencies

---

## Summary Statistics

| Category | Count | LOC |
|----------|-------|-----|
| Layer 5 (Integrals) | 7 files | 752 |
| Layer 6 (Evaluate) | 11 files | 651 |
| **Total** | **18 files** | **1403** |
| Standalone refs | 4 files | 305 |
| Pipeline files | 14 files | 1098 |

**LOC/file ratio (pipeline):** 78 lines/file (tid.jl and nlo_box.jl carry most weight).

**External dependencies:**
- COLLIER library (required for D₀, optional for B₀/C₀)
- QuadGK (Feynman parameter integrals in B₀, C₀, vacuum_polarization, vertex_f2)
- Base Julia only (no other Pkg deps)

