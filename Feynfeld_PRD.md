# Feynfeld.jl — Product Requirements Document

**Version:** 0.1-draft
**Date:** 2026-03-25
**Author:** Tobias J. Osborne / LUH
**Status:** Pre-development

---

## 1. Vision

Feynfeld.jl is a unified Julia package for symbolic and numerical quantum field theory computation. It replaces the fragmented Mathematica ecosystem of FeynRules, FeynArts, FeynCalc, FormCalc, and LoopTools with a single, end-to-end pipeline: from Lagrangian to cross-section in one `using Feynfeld`.

The fragmentation of the existing toolchain is an accident of Mathematica's limitations — namespace collisions between FeynArts and FeynCalc, Mathematica's inadequate numerical performance necessitating FORM and Fortran — not a reflection of logical boundaries in the physics. Feynfeld.jl dissolves these boundaries.

## 2. Motivation

### 2.1 Why this exists

No open-source, high-performance, end-to-end symbolic QFT package exists outside the Mathematica ecosystem. The existing tools are:

- **Fragmented**: Four packages (five counting FORM) that cannot coexist cleanly in one session.
- **Proprietary-dependent**: All require Mathematica licenses. FormCalc additionally requires FORM.
- **Slow**: Mathematica's symbolic engine is the bottleneck. FormCalc exists solely to offload work to FORM. LoopTools exists solely to offload numerics to Fortran.
- **Poorly tested**: Only FeynCalc has a proper test suite (MUnit). FeynArts and FormCalc have none.
- **Not composable**: Cannot be embedded in larger computational pipelines without Mathematica orchestration.

Julia eliminates all five problems simultaneously.

### 2.2 Immediate use case

The VLBAI experiment at Leibniz Universität Hannover searches for ultralight dark matter (ULDM) via atom interferometry. The theory pipeline from BSM Lagrangian → coupling constants → atomic sensitivity → interferometer phase response → exclusion limits currently requires manual stitching of Mathematica notebooks, Python scripts, and Fortran codes. Feynfeld.jl provides the first layer of this pipeline (Lagrangian → coupling constants) natively in Julia, enabling integration with downstream Julia codes for atomic physics, signal processing, and statistical analysis.

But ULDM is only the first testcase. The package is general-purpose.

### 2.3 Relationship to TensorGR.jl

TensorGR.jl (a Julia port of the xAct Mathematica suite for tensor calculus in general relativity) shares deep structural overlap with Feynfeld.jl:

- **Index contraction engine**: TensorGR.jl already implements abstract index notation, metric contraction, symmetry-aware canonicalization. Feynfeld.jl needs Lorentz index contraction over η^μν. These are the same algorithm over different metric signatures.
- **Tensor algebra**: Riemann tensor symmetries in TensorGR.jl and gamma-matrix algebra in Feynfeld.jl are both instances of term rewriting over indexed expressions with algebraic identities (Bianchi ↔ Clifford).
- **Code generation**: Both need to produce efficient numerical code from symbolic tensor expressions.

**Strategy**: Extract the shared index contraction and tensor algebra infrastructure from TensorGR.jl into a common foundation (or make TensorGR.jl a dependency). Feynfeld.jl's Lorentz layer is a specialisation of TensorGR.jl's abstract tensor layer to the Minkowski metric. This avoids reimplementing index gymnastics and gives both packages a shared, well-tested core.

Concretely:
- TensorGR.jl's `AbstractIndex`, `Metric`, `contract`, `canonicalize` → reused by Feynfeld.jl's Lorentz module.
- Feynfeld.jl adds: `DiracGamma`, `SUNMatrix`, `FourMomentum`, `SpinorChain` as new tensor-like types that compose with the shared index machinery.
- Gravitational coupling calculations (scalar field on curved background, graviton exchange) naturally live at the intersection: TensorGR.jl provides the geometry, Feynfeld.jl provides the QFT.

---

## 3. Architecture

### 3.1 One package, six layers

```
┌─────────────────────────────────────────────────┐
│                  Feynfeld.jl                     │
│                                                  │
│  ┌─────────────────────────────────────────────┐ │
│  │ Model                                       │ │
│  │ Lagrangian, fields, parameters, symmetries  │ │
│  └──────────────────┬──────────────────────────┘ │
│                     ▼                            │
│  ┌─────────────────────────────────────────────┐ │
│  │ Rules                                       │ │
│  │ Lagrangian → Feynman rules, vertex extract  │ │
│  └──────────────────┬──────────────────────────┘ │
│                     ▼                            │
│  ┌─────────────────────────────────────────────┐ │
│  │ Diagrams                                    │ │
│  │ Topology generation, field insertion,       │ │
│  │ amplitude construction                      │ │
│  └──────────────────┬──────────────────────────┘ │
│                     ▼                            │
│  ┌─────────────────────────────────────────────┐ │
│  │ Algebra  ← ← ← CORE, port first ← ← ← ←  │ │
│  │ Lorentz · Dirac · Colour · Tensor           │ │
│  │ Contraction, traces, simplify, PV reduction │ │
│  └──────────────────┬──────────────────────────┘ │
│                     ▼                            │
│  ┌─────────────────────────────────────────────┐ │
│  │ Integrals                                   │ │
│  │ PV scalar functions A₀ B₀ C₀ D₀            │ │
│  │ Symbolic + numerical evaluation             │ │
│  └──────────────────┬──────────────────────────┘ │
│                     ▼                            │
│  ┌─────────────────────────────────────────────┐ │
│  │ Evaluate                                    │ │
│  │ |M|², cross-sections, decay rates,          │ │
│  │ coupling RGE, codegen                       │ │
│  └─────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────┘
```

### 3.2 The Algebra layer is the type system

Everything above Algebra constructs expressions in its types. Everything below consumes them. The Algebra layer defines:

- `LorentzIndex`, `FourVector`, `FourMomentum`, `MetricTensor` — shared with / inherited from TensorGR.jl
- `DiracGamma`, `DiracChain`, `SpinorU`, `SpinorV`, `Slash` — Clifford algebra
- `SUNMatrix`, `SUNF`, `SUND`, `ColourDelta` — SU(N) colour algebra
- `PaVe`, `A0`, `B0`, `B1`, `C0`, `D0` — Passarino-Veltman integral symbols
- `Amplitude` — a sum of terms, each a product of the above

These are **not** Symbolics.jl expressions. They are a purpose-built algebraic type system with dispatch-based simplification. Mathematica's power for FeynCalc comes from pattern matching on symbolic expressions; Julia's equivalent is multiple dispatch on algebraic types.

### 3.3 Layer responsibilities

| Layer | ex-Package | Core operations |
|-------|-----------|-----------------|
| Model | FeynRules | Define fields, parameters, gauge groups; write Lagrangian as Julia expression |
| Rules | FeynRules | Functional differentiation of L w.r.t. fields → vertex functions |
| Diagrams | FeynArts | Enumerate topologies (graph theory), insert fields at Generic/Classes/Particles levels |
| Algebra | FeynCalc | Contract indices, Dirac traces, colour traces, tensor decomposition, PV reduction |
| Integrals | LoopTools | Numerical evaluation of scalar 1-loop integrals via 't Hooft-Veltman |
| Evaluate | FormCalc | Square amplitudes, sum/average over spins/colours, phase space integration |

FormCalc is **not ported**. Its reason for existence (Mathematica is slow, offload to FORM) does not apply. Its functionality (amplitude squaring, helicity method, Fortran codegen) is reimplemented natively in the Evaluate layer.

### 3.4 Dimensional regularisation

The algebra must work in both D and 4 dimensions throughout. Following FeynCalc's convention:
- D-dimensional objects carry explicit dimension tags
- The limit D → 4 is taken only at the end
- Support both conventional dimensional regularisation (CDR) and dimensional reduction (DRED) for SUSY compatibility

---

## 4. Build plan

### 4.1 Phase 0: Foundation (TensorGR.jl integration)

Extract or expose the abstract index contraction engine from TensorGR.jl. Define the `LorentzIndex` and `MetricTensor` types as specialisations. Validate: contracting g^μν g_νρ = δ^μ_ρ.

### 4.2 Phase 1: Algebra layer

Port FeynCalc's core functions, validated against FeynCalc's MUnit test suite.

Priority order (by dependency and test coverage):

1. **Lorentz algebra**: `Contract`, `MomentumExpand`, `MomentumCombine`, `Pair`, `ScalarProduct`
2. **Dirac algebra**: `DiracSimplify`, `DiracTrace`, `DiracOrder`, `SpinorChainTrick`
3. **Colour algebra**: `SUNSimplify`, `SUNTrace`, colour factor computation
4. **Tensor decomposition**: `TID` (tensor integral decomposition), `Tdec`
5. **PV reduction**: `PaVeReduce`, `PaVeOrder`, `PaVeAutoReduce`

Each function is ported, tested against the corresponding MUnit test, and only then integrated.

### 4.3 Phase 2: Integrals layer

Implement the scalar one-loop integrals A₀, B₀, B₁, B₀₀, C₀, D₀ numerically. Validate against LoopTools (call the Fortran library via `ccall` for comparison, then replace with native Julia).

### 4.4 Phase 3: Model + Rules layers

Implement the Lagrangian DSL and vertex extraction. Validate against FeynRules' SM model: extract all SM vertices and compare to published Feynman rules.

### 4.5 Phase 4: Diagrams layer

Implement topology generation and field insertion. Validate against FeynArts' example calculations (e+e- → ff̄ at tree level and one loop).

### 4.6 Phase 5: Evaluate layer

Amplitude squaring, spin/colour summation, phase-space integration, cross-section computation. Validate against textbook SM cross-sections (e+e- → μ+μ- at tree level as the first target).

### 4.7 Phase 6: ULDM application

Implement the scalar ULDM portal Lagrangian, extract coupling constants, compute sensitivity coefficients, connect to downstream VLBAI pipeline. This is the end-to-end integration test for the full stack.

---

## 5. Testing strategy

### 5.1 Ground truth hierarchy

1. **Published papers and textbooks**: Peskin & Schroeder cross-sections, Denner's one-loop formulae, PDG coupling constants. These are the ultimate ground truth.
2. **FeynCalc MUnit tests**: The primary porting oracle. Each MUnit test is translated to a Julia `@test` and must pass before the corresponding function is considered ported.
3. **Cross-validation against reference implementations**: For numerical results, compare against LoopTools (Fortran), COLLIER, or Package-X. At least two independent checks per numerical function.
4. **Physics invariants**: Ward identities, gauge invariance of physical observables, unitarity cuts. These are structural tests that no specific implementation can fake.

**Rule 6 applies absolutely**: physics is ground truth, not pinned numbers. If a FeynCalc MUnit test and a textbook disagree, the textbook wins and the test is suspect. Investigate before proceeding.

### 5.2 Test infrastructure

- Every exported function has at least one `@test`.
- Tests are organized to mirror the layer structure: `test/algebra/`, `test/integrals/`, `test/model/`, etc.
- A `test/references/` directory contains local copies of published formulae (from papers, not from LLM memory) used as ground truth.
- CI runs the full test suite. No PR merges without green tests.

### 5.3 MUnit translation protocol

For each FeynCalc MUnit test file:
1. Obtain a local copy of the MUnit `.mt` file from the FeynCalc repository.
2. Translate each `MUnit`Test` to a Julia `@test`, preserving the mathematical content exactly.
3. Where the MUnit test checks symbolic equality, the Julia test checks algebraic equivalence (expressions may be in different canonical forms).
4. Where the MUnit test checks numerical equality, the Julia test checks to the same tolerance.
5. Document the source MUnit file and test ID in a comment above each `@test`.

---

## 6. Development rules

These are non-negotiable. They are learnt through suffering.

### TOBIAS'S RULES — FOLLOW TO THE LETTER

1. **SKEPTICISM**: All subagent work, handoffs — verify everything twice.
2. **DEEP BUGS**: Deep, complex, interlocked. Do not underestimate.
3. **NO BANDAIDS**: Best-practices full solutions only.
4. **WORKFLOW**: 3 subagents before any core code change (research local copy of port source + 2 solutions).
5. **REVIEW**: Rigorous reviewer agent after every core change. No exceptions.
6. **GROUND TRUTH**: Physics is ground truth, not pinned numbers. Tests may be suspect. Local copies of published papers, or reference implementations are the ONLY truth. Everything else is hallucination.
7. **TESTING**: Targeted only, or full suite in background.
8. **REPEAT RULES**: Repeat occasionally to maintain focus.
9. **DO NOT UNDERESTIMATE**: This is deeply nontrivial.
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Run agents sequentially only.

### Additional development constraints

- **LOC limit**: No source file exceeds ~200 lines. Split aggressively.
- **Handoff documentation**: Every work session ends with a handoff document describing what was done, what works, what doesn't, and what to do next.
- **Beads issue tracker**: Use the Beads system for tracking issues. Each bead is a self-contained unit of work with clear entry/exit criteria.
- **Anti-panic protocol**: When a deep bug is encountered, stop. Document. Research. Do not thrash. Three subagent investigations before any fix attempt.

---

## 7. Design principles

### 7.1 Types over patterns

Mathematica FeynCalc uses rule-based pattern matching on symbolic trees. Julia Feynfeld uses multiple dispatch on algebraic types. The translation principle:

- Mathematica `f[x_?SomeTest] := ...` → Julia `f(x::SomeType) = ...`
- Mathematica `expr /. {rule1, rule2}` → Julia method dispatch or explicit `simplify` passes
- Mathematica `Hold`, `HoldForm` → Julia expression types with lazy evaluation

### 7.2 Immutable expressions, functional transformations

All algebraic expressions are immutable. Simplification produces new expressions, never mutates. This enables:
- Safe parallel evaluation (future)
- Reliable equality testing
- Clear debugging (print any intermediate state)

### 7.3 Canonical forms

Every expression type has a unique canonical form. Two mathematically equal expressions must produce identical canonical forms. This is tested exhaustively. The canonical form is the normal-ordered, fully-contracted, colour-decomposed representation.

### 7.4 Performance is not premature

Julia is chosen partly for performance. The algebra layer should be fast enough that FormCalc/FORM is unnecessary. Profile early. If Dirac trace computation is slower than FeynCalc+Mathematica, something is wrong.

### 7.5 Interop

- **TensorGR.jl**: Shared index algebra foundation. GR + QFT calculations in one session.
- **Symbolics.jl**: Optional backend for symbolic parameter manipulation (masses, couplings as symbolic variables). Not used for the core algebra — too general, too slow for index-heavy computation.
- **LoopTools via ccall**: During development, call the Fortran LoopTools library for numerical validation. Replace with native Julia implementation once validated.

---

## 8. Scope boundaries

### 8.1 In scope (v1.0)

- Full SM at tree level and one loop
- Arbitrary BSM models via Lagrangian input
- Scalar, fermion, vector field types
- SU(N) gauge groups
- Dimensional regularisation (CDR and DRED)
- One-loop PV reduction and numerical evaluation
- Cross-section and decay rate computation for 2→2 processes
- UFO-compatible model file import (read existing FeynRules models)

### 8.2 Out of scope (v1.0, future work)

- Multi-loop integrals (IBP reduction, sector decomposition)
- Real radiation and IR subtraction (FKS, Catani-Seymour)
- Parton distribution functions and hadron-level cross-sections
- SUSY-specific features beyond DRED
- Spin-3/2, spin-2 fields
- Automated renormalisation (except for specific models)
- Monte Carlo event generation

### 8.3 Out of scope (permanently)

- Reimplementing Mathematica's general-purpose CAS
- GUI or interactive topology editor
- Backward compatibility with FeynCalc's Mathematica API

---

## 9. Success criteria

### 9.1 Minimum viable

Feynfeld.jl can compute the tree-level cross-section for e+e- → μ+μ- in QED starting from the QED Lagrangian, and the result agrees with Peskin & Schroeder equation (5.10) to machine precision.

### 9.2 ULDM milestone

Feynfeld.jl can derive the ULDM-electron and ULDM-photon coupling constants (d_{m_e}, d_e) from the scalar portal Lagrangian L_φ = φ√(4πG_N) [d_e/(4e²) F_μν F^μν − d_{m_e} m_e ψ̄_e ψ_e], compute the one-loop radiative corrections, and produce numerical coupling constants suitable for input to the VLBAI signal pipeline.

### 9.3 Community milestone

Feynfeld.jl can reproduce the sensitivity projections from Badurina et al. (2109.10965) for AION-10/100/km scalar ULDM searches, matching published figures to within plotting accuracy.

### 9.4 Full SM milestone

Feynfeld.jl can compute all 2→2 SM processes at one loop, matching published results from Denner et al. and/or the FormCalc HEP Process Repository.

---

## 10. Risks

| Risk | Severity | Mitigation |
|------|----------|------------|
| Dirac algebra canonicalization is subtler than expected | High | Rule 9. Port FeynCalc's approach exactly first, optimise later. |
| FeynCalc MUnit tests contain bugs | Medium | Rule 6. Cross-validate against textbooks. |
| TensorGR.jl integration introduces coupling | Medium | Clean interface boundary. Depend on abstract types only. |
| Colour algebra for exceptional groups (E₆, E₈) | Low | Not needed for SM or ULDM. Defer. |
| Performance regression vs Mathematica+FORM | Medium | Profile continuously. Julia should win; if not, investigate. |
| Scope creep into multi-loop | High | Hard scope boundary at v1.0. One loop only. |

---

## 11. References

### Primary sources for validation

- M. E. Peskin, D. V. Schroeder, *An Introduction to Quantum Field Theory* (1995)
- A. Denner, "Techniques for the Calculation of Electroweak Radiative Corrections at the One-Loop Level and Results for W-physics at LEP2", Fortschr. Phys. 41 (1993) 307
- V. Shtabovenko, R. Mertig, F. Orellana, "FeynCalc 10", arXiv:2312.14089
- T. Hahn, "Generating Feynman Diagrams and Amplitudes with FeynArts 3", hep-ph/0012260
- T. Hahn, M. Pérez-Victoria, "FormCalc", Comput. Phys. Commun. 118 (1999) 153
- A. Alloul et al., "FeynRules 2.0", Comput. Phys. Commun. 185 (2014) 2250
- G. 't Hooft, M. Veltman, "Scalar One-Loop Integrals", Nucl. Phys. B153 (1979) 365
- G. Passarino, M. Veltman, "One-Loop Corrections for e+e- Annihilation into μ+μ- in the Weinberg Model", Nucl. Phys. B160 (1979) 151

### ULDM-specific

- L. Badurina et al., "Refined ultralight scalar dark matter searches with compact atom gradiometers", Phys. Rev. D 105 (2022) 023006, arXiv:2109.10965
- L. Badurina et al., "Prospective Sensitivities of Atom Interferometers to GWs and ULDM", arXiv:2108.02468
- L. Badurina et al., "AION: An Atom Interferometer Observatory and Network", arXiv:1911.11755
- Centers et al., "Stochastic fluctuations of bosonic dark matter", Nat. Commun. 12 (2021) 7321

---

## Appendix A: Name

*Feynfeld* = Feynman + Feld (German: field). Consistent with the -feld naming convention: Abstractfeld.jl, Tensorfeld, Convexfeld, Alethfeld.
