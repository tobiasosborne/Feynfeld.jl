# Spiral 9 Implementation Plan: Architecture Consolidation

**Date:** 2026-03-29
**Status:** Plan (research complete, design validated by 2 agents)
**Goal:** Every existing process flows through the full pipeline. No more hand-built amplitudes.
**Scope:** ~200 LOC new code across 3 files. NOT a FeynArts port.

---

## The Problem

Layer 4 (Algebra) is excellent. But Layers 1-3 are scaffolding:

```
CURRENT STATE:
  e+e-→μ+μ-    → Model → Rules → Diagrams → Algebra → σ     ✓ full pipeline
  Compton       → [hand-built in test file]  → Algebra → σ     ✗ bypasses L1-3
  Bhabha        → [hand-built in test file]  → Algebra → σ     ✗ bypasses L1-3
  qq̄→gg        → [hand-built in test file]  → Algebra → σ     ✗ bypasses L1-3
  e+e-→W+W-    → [standalone formula]        → σ               ✗ bypasses L1-5
```

The pipeline is validated for 1 of 7 processes. This ossifies if we don't fix it now.

---

## The Solution: 3 files, ~200 LOC

**Key insight (both research agents converge):** Tree-level 2→2 has exactly
3 channel topologies (s, t, u). Diagram generation is FILTERING, not graph
enumeration. No external library needed.

### File 1: `channels.jl` (~70 LOC)

```julia
# A tree diagram is fully specified by 4 legs + 1 internal line + 2 vertices.
struct TreeDiagram22
    channel::Symbol              # :s, :t, :u
    exchanged::Symbol            # field name of internal propagator
    vertex_L::NTuple{3, Symbol}  # fields at left vertex
    vertex_R::NTuple{3, Symbol}  # fields at right vertex
    propagator_mom::Momentum     # or MomentumSum
end

# Enumerate valid diagrams by checking vertex existence in the model.
# Complexity: O(3 channels × N_fields) ≈ 60 checks for SM.
function tree_diagrams_22(model, rules, in_fields, out_fields)
    diagrams = TreeDiagram22[]
    for internal in model_fields(model)
        # s-channel: (in1,in2) → internal → (out1,out2)
        # t-channel: (in1,out1) via internal, (in2,out2) via internal
        # u-channel: (in1,out2) via internal, (in2,out1) via internal
        # Check vertex existence in rules for each configuration
    end
    diagrams
end
```

### File 2: `amplitude.jl` (~80 LOC)

Generic amplitude builder — works for ANY TreeDiagram22.
Dispatch happens at the vertex/propagator level (already exists), not at channel level.

```julia
function build_amplitude(d::TreeDiagram22, rules, momenta, spinors)
    # Look up vertex Lorentz structures via existing dispatch:
    #   vertex_structure(::Fermion, ::Fermion, ::Boson) → gamma chain
    #   vertex_structure(::Boson, ::Boson, ::Boson) → triple gauge tensor
    # Look up propagator via existing dispatch:
    #   propagator_num(::Fermion, p, m) → (p̸ + m)
    #   propagator_num(::Boson, p, M) → -g^μν or -g^μν + k^μk^ν/M²
    # Wire together: spinor × vertex_L × propagator × vertex_R × spinor
end
```

Two cases determined by exchanged field species:
- **Boson exchange** (ee→μμ s, Bhabha s/t): 2 separate DiracChains contracted by metric
- **Fermion exchange** (Compton s/u): 1 DiracChain with propagator numerator inside

### File 3: `ew_model.jl` (~60 LOC)

Explicit SM vertex table — NOT a Lagrangian parser.

```julia
struct EWModel <: AbstractModel
    # reuses ew_parameters.jl constants
end

# Vertex table as dispatch. The physics lives HERE.
# ffV: standard vector coupling
vertex_structure(::EWModel, ::Fermion, ::Fermion, ::Photon, mu) = # -ieQγ^μ
# ffVA: chiral coupling (needs gamma5)
vertex_structure(::EWModel, ::Fermion, ::Fermion, ::ZBoson, mu) = # -ie/(s_Wc_W)(g_V-g_Aγ5)γ^μ
# VVV: triple gauge coupling (momentum-dependent)
vertex_structure(::EWModel, ::WPlus, ::WMinus, ::Photon, mu, nu, rho) = # standard WWγ
```

---

## Dependency Chain

```
Spiral 8 (current):
  [bug fixes] ─┬─ colour_trace i²
               ├─ B₀ imaginary
               ├─ DiracExpr.+ growth
               └─ Eps contraction ──→ gamma5 traces
                                        │
Spiral 9:                               │
  Channel type (feynfeld-e8d)            │
       │                                │
  build_amplitude (feynfeld-czy)         │
       │                                │
  Pipeline tests (feynfeld-07g)     EW model (feynfeld-cdf)
       │                                │
  Mark standalone refs (feynfeld-cur)    │
                                        │
Spiral 10:                              │
  D₀ evaluation ──→ box diagrams        │
                                        │
Spiral 11+:                             │
  Lagrangian DSL ──→ BSM models ────────┘
```

---

## Implementation Phasing (from research)

**Phase A: QED boson exchange** (~60 LOC, covers ee→μμ + Bhabha s-channel)
- TreeDiagram22 struct
- tree_diagrams_22 for boson-exchange channels
- build_amplitude for boson exchange (2 separate DiracChains)
- Pipeline test: verify ee→μμ matches existing test_vertical.jl result

**Phase B: QED fermion exchange** (~40 LOC, covers Compton + Bhabha t-channel)
- build_amplitude for fermion exchange (1 DiracChain with p̸+m)
- Interference term handling (Bhabha: 8-gamma single trace for cross terms)
- Pipeline test: verify Compton and Bhabha match existing results

**Phase C: QCD** (~40 LOC, covers qq̄→gg)
- QCDModel with qqg vertex + ggg triple vertex
- Axial gauge polarization sum integration
- Analytical colour factors (already implemented)
- Pipeline test: verify qq̄→gg matches existing result

**Phase D: EW** (~60 LOC, covers ee→W+W-, the hardest)
- EWModel with eeγ, eeZ, eνW, WWγ, WWZ vertices
- Requires gamma5 traces (blocked on feynfeld-d2g)
- Massive vector propagator and polarization sum (already implemented)
- Pipeline test: verify ee→WW matches existing Grozin formula

---

## What NOT to build

1. **No graph library.** Tree-level 2→2 = {s, t, u}. Enumerate, don't generate.
2. **No Lagrangian parser.** SM vertices are known. Hard-code them as dispatch methods.
3. **No .gen/.mod files.** Julia types ARE the model definition.
4. **No symbolic manipulation.** Vertex factors are concrete DiracChain builders, not symbolic rules.
5. **No abstraction for hypothetical futures.** Build for SM tree+1-loop. BSM comes later.

---

## Success Criteria

Spiral 9 is DONE when:

```julia
# This works for ALL 7 processes, through the full pipeline:
model = ew_model()
channels = tree_channels(model, [:e, :ebar], [:mu, :mubar])
amps = [build_amplitude(ch, feynman_rules(model), momenta) for ch in channels]
m_squared = spin_sum_amplitude_squared(amps...)
sigma = evaluate_cross_section(m_squared, sqrt_s)

# And the result matches the existing hand-built / analytical values
@test sigma ≈ known_sigma rtol=1e-6
```

No more bypassing the pipeline. Every process goes through Model → Channels → Amplitude → Algebra → σ.

---

## Beads Issues

| ID | Title | Priority | Blocked by |
|----|-------|----------|------------|
| feynfeld-e8d | Channel type + tree_channels() | P1 | nothing |
| feynfeld-czy | build_amplitude() dispatch | P1 | feynfeld-e8d |
| feynfeld-cdf | EWModel with SM vertices | P1 | feynfeld-d2g (gamma5) |
| feynfeld-2yw | QCDModel with gluon vertices | P2 | nothing |
| feynfeld-07g | Pipeline tests for all processes | P1 | feynfeld-czy |
| feynfeld-cur | Mark standalone recipes as reference | P2 | feynfeld-07g |
| feynfeld-ntj.8 | SM model file | P2 | feynfeld-cdf |
| feynfeld-ntj.9 | Validate SM vertices | P2 | feynfeld-ntj.8 |
