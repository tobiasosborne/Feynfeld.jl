# Vertical Tracer Bullet Plan: e+e- → μ+μ- through all 6 layers

## Process: e+e- → μ+μ- (QED, tree + 1-loop)

### Why this process?
- Algebra layer already validated (P&S Eq 5.10): 8(t² + u²)
- Self-energy DiracTrick validated in D dimensions
- Known textbook answer at every intermediate step
- Tree exercises layers 1-2-3-4-6; 1-loop adds layer 5
- Single unified test, no wasted work

---

## LAYER 1: Model (~60 LOC)

### Input: Nothing (user-facing entry point)
### Output: `QEDModel` struct containing fields, couplings, parameters

### Design:

```julia
# Fields as parametric types
abstract type FieldType end
struct FermionField <: FieldType
    name::Symbol
    mass::Symbol        # symbolic mass parameter
    charge::Rational{Int}  # electric charge in units of e
end
struct VectorField <: FieldType
    name::Symbol
    mass::Symbol
end

# A Model is a collection of fields + coupling constants
struct QFTModel
    name::Symbol
    fields::Vector{FieldType}
    params::Dict{Symbol, Any}  # :m_e => 0.511 MeV, :e => coupling, etc.
    gauge_group::Symbol         # :U1, :SU3, etc.
end

# QED convenience constructor
function qed_model(; m_e=:m_e, m_mu=:m_mu)
    QFTModel(:QED,
        [FermionField(:e, m_e, -1//1),
         FermionField(:mu, m_mu, -1//1),
         VectorField(:gamma, :zero)],
        Dict(:m_e => m_e, :m_mu => m_mu, :e => :e, :alpha => :(e^2/(4*pi))),
        :U1
    )
end
```

### Validation:
- `model = qed_model()` constructs without error
- `model.fields` has 3 entries (e, mu, gamma)
- `model.gauge_group == :U1`

### File: `src/v2/model.jl`

---

## LAYER 2: Rules (~100 LOC)

### Input: `QFTModel`
### Output: `FeynmanRules` struct containing vertices + propagators

### Design:

```julia
# A vertex is: (field_names...) → DiracExpr (the coupling structure)
struct VertexRule
    fields::Vector{Symbol}           # e.g., [:e, :e, :gamma]
    lorentz_structure::Function       # (indices...) → DiracExpr
    coupling::Any                     # :(-im*e) or similar
end

# A propagator is: field → (momentum → expression)
struct PropagatorRule
    field::Symbol
    structure::Function  # (p, mass) → DiracExpr or AlgSum
end

struct FeynmanRules
    vertices::Vector{VertexRule}
    propagators::Vector{PropagatorRule}
end

# Extract QED Feynman rules
function feynman_rules(model::QFTModel)
    @assert model.gauge_group == :U1 "Only U(1) supported"

    # QED vertex: fermion-fermion-photon
    # Rule: -ie γ^μ  (μ is the photon Lorentz index)
    fermion_fields = [f for f in model.fields if f isa FermionField]
    vertices = VertexRule[]
    for f in fermion_fields
        push!(vertices, VertexRule(
            [f.name, f.name, :gamma],
            (mu,) -> DiracExpr(DiracChain([DiracGamma(LISlot(mu))])),
            :(-im*e)  # coupling (i from Feynman rule, e from charge)
        ))
    end

    # Propagators
    propagators = [
        # Fermion: i(p-slash + m) / (p² - m²)
        PropagatorRule(:fermion, (p, mass) -> ...),
        # Photon: -i g^{μν} / p²  (Feynman gauge)
        PropagatorRule(:photon, (p, mu, nu) -> ...),
    ]

    FeynmanRules(vertices, propagators)
end
```

### Validation:
- `rules = feynman_rules(qed_model())`
- `length(rules.vertices) == 2` (e-e-γ and μ-μ-γ)
- Vertex structure matches P&S Table of Feynman rules

### File: `src/v2/rules.jl`

---

## LAYER 3: Diagrams (~120 LOC)

### Input: `FeynmanRules` + process specification (in/out particles)
### Output: Vector of `FeynmanDiagram` structs

### Design:

For the tracer bullet, we hard-code the single tree-level diagram.
Full topology generation is deferred.

```julia
struct ExternalLeg
    field::Symbol
    momentum::Momentum
    incoming::Bool
end

struct InternalLine
    field::Symbol
    momentum::Momentum  # or MomentumSum
end

struct FeynmanDiagram
    external::Vector{ExternalLeg}
    internal::Vector{InternalLine}
    vertices::Vector{Symbol}           # vertex labels
    amplitude::DiracExpr               # the unsimplified amplitude
end

# Hard-coded tree-level e+e- → μ+μ- diagram generator
function tree_diagrams_ee_mumu(rules::FeynmanRules)
    # Single s-channel diagram: e+e- → γ* → μ+μ-
    #
    #   e-(p1) ──→──┐     ┌──→── μ-(k1)
    #               ├─γ*──┤
    #   e+(p2) ──→──┘     └──→── μ+(k2)
    #
    # Amplitude = [v̄(p2) (-ieγ^μ) u(p1)] × (-ig_{μν}/s) × [ū(k1) (-ieγ^ν) v(k2)]
    #           = (ie)² / s × [v̄(p2) γ^μ u(p1)] [ū(k1) γ_μ v(k2)]

    p1 = Momentum(:p1); p2 = Momentum(:p2)
    k1 = Momentum(:k1); k2 = Momentum(:k2)
    mu = LorentzIndex(:mu)

    external = [
        ExternalLeg(:e, p1, true),   # incoming e-
        ExternalLeg(:e, p2, true),   # incoming e+ (antiparticle)
        ExternalLeg(:mu, k1, false), # outgoing μ-
        ExternalLeg(:mu, k2, false), # outgoing μ+ (antiparticle)
    ]

    # Build amplitude chains
    chain1 = dot(vbar(p2), GA(:mu), u(p1))     # electron current
    chain2 = dot(ubar(k1), GA(:mu), v(k2))     # muon current

    [FeynmanDiagram(external, [...], [:eeγ, :μμγ], ...)]
end
```

### Validation:
- `diagrams = tree_diagrams_ee_mumu(rules)`
- `length(diagrams) == 1` (single s-channel)
- Diagram has 4 external legs, correct momenta
- Amplitude structure matches hand-written tracer bullet

### Files: `src/v2/diagrams.jl`

### For 1-loop (Stage B), add:
- Vertex correction diagram (triangle)
- Self-energy correction diagrams (2)
- Vacuum polarization diagram (1)
- Each with loop momentum `k` and FeynAmpDenominator

---

## LAYER 4: Algebra (DONE ✓)

### Input: `FeynmanDiagram.amplitude`
### Output: Simplified `AlgSum` (all indices contracted, scalar result)

### Pipeline:
```julia
function simplify_amplitude(diag::FeynmanDiagram; ctx::SPContext)
    # 1. Spin sum (completeness relations)
    traced = spin_sum_and_trace(diag)    # → AlgSum

    # 2. Contract Lorentz indices
    contracted = contract(traced)         # → AlgSum

    # 3. Expand scalar products
    expanded = expand_scalar_product(contracted)  # → AlgSum

    # 4. Evaluate SP in terms of Mandelstam variables
    evaluate_sp(expanded; ctx=ctx)        # → AlgSum (pure number)
end
```

### Validation: P&S Eq (5.10): `|M|² = 8(t² + u²)`
### Already validated in `test/v2/test_ee_mumu_x.jl`

---

## LAYER 5: Integrals (~200 LOC, Stage B only)

### Input: Loop amplitude with `FeynAmpDenominator`
### Output: Numerical values of scalar 1-loop integrals

### What's needed:

**5a. Tensor decomposition** (v1 has tdec for rank 0-2, port to v2):
```julia
# Integral of k^μ / [k²((k-p)²-m²)] = p^μ B1(p²; 0, m²)
tdec(rank, ext_momenta, masses) → PaVe coefficients
```

**5b. PaVe reduction** (v1 has B-function reduction):
```julia
# B1(p²; m1², m2²) = [A0(m2²) - A0(m1²) - (p²+m1²-m2²)B0] / (2p²)
pave_reduce(B1(...)) → expression in A0, B0
```

**5c. Numerical evaluation** (NEW — the big piece):
```julia
# A0(m²) = m²(1 - log(m²/μ²))  [in dim-reg, MS-bar]
# B0(p²; m1², m2²) = ∫₀¹ dx log(...)  [Feynman parameter integral]
# Need: Li₂(x) = -∫₀ˣ log(1-t)/t dt  (dilogarithm)

function A0_numerical(m²; μ²=1.0)::ComplexF64
    m² * (1 - log(m²/μ²))  # + O(ε) pole handled separately
end

function B0_numerical(p², m1², m2²; μ²=1.0)::ComplexF64
    # Feynman parameter integration or 't Hooft-Veltman formula
    ...
end
```

### Validation:
- `A0(m²)` matches LoopTools `A0i(aa0, m²)`
- `B0(s, 0, m²)` matches LoopTools `B0i(bb0, s, 0, m²)`
- For the self-energy: `B0(p², 0, m_e²)` at benchmark kinematics

### Files:
- `src/v2/pave_types.jl` — PaVe{N} type (port from v1)
- `src/v2/pave_reduce.jl` — B-function reduction
- `src/v2/pave_numerical.jl` — A0, B0, Li2 numerical evaluation

### Reference: `/refs/LoopTools/` Fortran source, 't Hooft-Veltman (1979)

---

## LAYER 6: Evaluate (~150 LOC)

### Input: Simplified `|M|²` from Layer 4, optionally with loop corrections
### Output: Differential/total cross-section

### Design:

```julia
# Phase space for 2 → 2 massless particles
struct TwoBodyPhaseSpace
    s::Number   # CM energy squared
    cosθ::Number  # scattering angle
end

# Mandelstam variables from phase space
function mandelstam(ps::TwoBodyPhaseSpace)
    s = ps.s
    t = -s/2 * (1 - ps.cosθ)
    u = -s/2 * (1 + ps.cosθ)
    (s=s, t=t, u=u)
end

# Differential cross-section: dσ/dΩ = |M|² / (64π²s)
function dsigma_domega(M_squared::Number, s::Number)
    M_squared / (64 * π^2 * s)
end

# Total cross-section: integrate dσ/dΩ over solid angle
function sigma_total(M_squared_func, s::Number)
    # M_squared_func(cosθ) → |M|²
    # σ = ∫₋₁¹ dcosθ ∫₀²π dφ × dσ/dΩ
    # For azimuthal symmetry: σ = 2π ∫₋₁¹ dcosθ × dσ/dΩ
    # Gauss-Legendre quadrature or analytical for tree-level
    ...
end

# Coupling: include e⁴/s² prefactor
# Full: dσ/dΩ = α²/(4s) × (1 + cos²θ)  [P&S Eq 5.11]
```

### Validation (tree-level):
- `dσ/dΩ(θ=π/2) = α²/(4s)` (P&S)
- `σ_total = 4πα²/(3s)` (P&S Eq 5.12, massless limit)
- At √s = 91.2 GeV (Z pole): σ ≈ 87 nb (known experimental value gives sanity check)

### Validation (1-loop, Stage B):
- O(α) correction: `δσ/σ = (3α)/(4π) × [π²/3 - 1/2]` (Schwinger correction)
- This is one of the most famous results in QED

### Files:
- `src/v2/cross_section.jl` — phase space + dσ/dΩ
- `src/v2/polarization_sum.jl` — photon polarization sums (for Compton later)

---

## THE TEST FILE

```julia
# test/v2/test_vertical.jl — the full vertical tracer bullet
#
# Exercises ALL layers of Feynfeld for e+e- → μ+μ-:
#   Model → Rules → Diagrams → Algebra → Evaluate
#
# Stage A: tree-level → P&S (5.10) and (5.12)
# Stage B: 1-loop → Schwinger correction (future)

@testset "Vertical tracer bullet: e+e- → μ+μ-" begin

    @testset "Layer 1: Model" begin
        model = qed_model()
        @test length(model.fields) == 3
        @test model.gauge_group == :U1
    end

    @testset "Layer 2: Rules" begin
        rules = feynman_rules(qed_model())
        @test length(rules.vertices) >= 2   # e-e-γ and μ-μ-γ
    end

    @testset "Layer 3: Diagrams" begin
        rules = feynman_rules(qed_model())
        diagrams = tree_diagrams_ee_mumu(rules)
        @test length(diagrams) == 1  # single s-channel
    end

    @testset "Layer 4: Algebra (already validated)" begin
        # Reuse existing test_ee_mumu_x.jl logic
        # Verify: |M|² = 8(t² + u²) at benchmark kinematics
    end

    @testset "Layer 6: Cross-section (tree)" begin
        # dσ/dΩ = α²/(4s) × (1 + cos²θ)
        # σ_total = 4πα²/(3s)
        α = 1/137.036
        s = (91.2)^2  # GeV², Z-pole energy
        σ = 4π * α^2 / (3s)  # in GeV^{-2}
        # Convert: 1 GeV^{-2} = 0.3894 mb = 389.4 μb = 389400 nb
        σ_nb = σ * 389400  # nanobarns
        @test σ_nb ≈ 1.45  atol=0.1  # ~1.45 nb at Z-pole (QED only, no Z)
    end

    # @testset "Layer 5: Integrals (Stage B)" begin
    #     # Deferred to Stage B
    #     # B0(s, 0, m_e²) at benchmark kinematics
    # end

    # @testset "Stage B: 1-loop correction" begin
    #     # Schwinger correction: δσ/σ = (3α)/(4π)(π²/3 - 1/2) ≈ 0.17%
    # end
end
```

---

## IMPLEMENTATION ORDER

### Stage A (Tree-level, exercises Layers 1-2-3-4-6):

| Step | Layer | File | LOC | What |
|------|-------|------|-----|------|
| 1 | Model | `src/v2/model.jl` | ~60 | QFTModel, FermionField, VectorField, qed_model() |
| 2 | Rules | `src/v2/rules.jl` | ~100 | FeynmanRules, vertices, propagators |
| 3 | Diagrams | `src/v2/diagrams.jl` | ~120 | FeynmanDiagram, tree_diagrams_ee_mumu() |
| 4 | Evaluate | `src/v2/cross_section.jl` | ~80 | Phase space, dσ/dΩ, σ_total |
| 5 | Test | `test/v2/test_vertical.jl` | ~150 | Full pipeline test |
| | | **Total** | **~510** | |

### Stage B (1-loop, adds Layer 5):

| Step | Layer | File | LOC | What |
|------|-------|------|-----|------|
| 6 | Integrals | `src/v2/pave_types.jl` | ~60 | PaVe{N} type system |
| 7 | Integrals | `src/v2/pave_reduce.jl` | ~80 | B-function reduction |
| 8 | Integrals | `src/v2/pave_numerical.jl` | ~120 | Li₂, A₀, B₀ numerical |
| 9 | Diagrams | `src/v2/diagrams.jl` ext | ~80 | Loop diagram generation |
| 10 | Test | `test/v2/test_vertical.jl` ext | ~100 | Stage B tests |
| | | **Total** | **~440** | |

### Grand total: ~950 new LOC for all 6 layers

---

## OPEN DESIGN QUESTIONS (resolve before coding)

1. **Model DSL**: `qed_model()` convenience vs `@model QED begin ... end` macro?
   - Recommendation: Start with plain struct constructors. Add macros later.

2. **Rules extraction**: Tobias requires second-quantization operator algebra,
   NOT functional differentiation. For QED tree-level this is equivalent.
   What matters: the INTERFACE (FeynmanRules struct) is the same either way.
   - Recommendation: Hard-code QED rules for tracer bullet, implement
     second-quantization extraction as a separate step later.

3. **Diagram generation**: Full recursive graph enumeration vs hard-coded?
   - Recommendation: Hard-code for tracer bullet. Layer 3's VALUE is the
     interface (FeynmanDiagram struct), not the algorithm.

4. **Numerical integrals**: Native Julia or FFI to LoopTools Fortran?
   - Recommendation: Native Julia (pure, portable, testable).
     Li₂ is ~20 LOC. A0/B0 are ~50 LOC each. No Fortran needed.

5. **Phase space integration**: Analytical (for 2→2 tree) or numerical?
   - Recommendation: Analytical for tree (exact). Numerical for loop corrections.
