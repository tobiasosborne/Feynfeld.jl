# FeynfeldX Design Choices — Lessons from the v2 Experiment

## READ THIS FIRST. This is the onboarding doc for the v2 experimental rebuild.

## Why v2 exists

v1 mirrored FeynCalc's Mathematica patterns: tagged unions, `Expr`-tree coefficients,
polymorphic return types. This imported Mathematica's weaknesses into Julia.
v2 asks "how would a Julia expert design this?" and goes hard on multiple dispatch,
parametric types, and type stability. The tracer bullet (e+e- → mu+mu-) took 3 sessions
in v1 and 30 minutes in v2. The full vertical test (Model → Rules → Diagrams → Algebra
→ Evaluate) was built and passing in a single session.

---

## Core Principles (VALIDATED)

1. **Parametric types, not tagged unions.** If you write `x.field isa SomeType`, you're
   doing it wrong. Make the type a parameter and dispatch.
2. **Uniform return types.** Scalar algebra returns `AlgSum`. Matrix-valued Dirac
   expressions return `DiracExpr`. No `Int|Pair|Eps|Expr|Vector` polymorphism.
3. **Proper coefficient algebra.** Coefficients are `DimPoly` (polynomial in D) or
   `Rational{Int}`. Never raw `Expr`. Never `Any`.
4. **Dict-based expression storage.** `AlgSum` is `Dict{FactorKey, Coeff}`. Like-term
   collection is O(1).
5. **Dispatch is the control flow.** Physics spaces (Lorentz, Dirac, colour) are
   encoded in types. `contract`, `trace`, `simplify` dispatch on them.
6. **Construction is cheap, simplification is explicit.** Constructors canonicalise
   (sort, normalise signs) but never simplify.
7. **Research Julia idioms before every new layer.** Don't default to OOP struct
   hierarchies. Check how Julia ecosystem packages solve the same problem.

---

## Lessons Learned (Session 4 — experiment session)

### What worked well

| Pattern | Where used | Why it worked |
|---------|-----------|---------------|
| **Parametric types** | `Pair{A,B}`, `DiracGamma{S}`, `Spinor{K}`, `Field{Species}` | Eliminated all `isa` cascades. Compiler does the dispatch. |
| **DimPoly coefficients** | `coeff.jl` | Killed the entire `_mul_coeff`/`_flatten_product`/`_to_algsum` machinery from v1. One type, clean arithmetic. |
| **Dict-based AlgSum** | `expr.jl` | Like-term collection is free. No `_collect_terms` scan. |
| **Uniform return types** | `dirac_trace → AlgSum`, `colour_trace → AlgSum` | No adaptation layers. Functions compose directly. |
| **Abstract interface** | `AbstractModel`, `model_fields()`, `gauge_groups()` | Multiple model backends possible without changing downstream code. |
| **Gauge groups as types** | `U1`, `SU{N}` | Casimirs, dimensions as compile-time dispatch. |
| **Field species as type param** | `Field{Fermion}`, `Field{Boson}` | Propagator dispatch: `propagator_num(::Fermion, ...)` vs `propagator_num(::Boson, ...)` |
| **Callable structs** | `FeynmanRules((:e,:e,:gamma))` | Natural pipeline. Each layer carries its state. |
| **ScopedValues** | `SPContext` via `CURRENT_SP` | Clean implicit context, explicit override for testing. |
| **Expand-before-trace** | `MomSumSlot` handling | Don't complicate core algorithms. Pre-process inputs. |

### What caused friction (COCKROACHES)

| Friction point | What broke | Root cause | Fix/lesson |
|----------------|-----------|------------|------------|
| **`Base.Pair` collision** | `sp_context(...)` couldn't accept `=>` syntax | Our `Pair` shadows `Base.Pair` | **Don't export `Pair`**. Users use `SP/FV/MT`. Internal code uses `Base.Pair` explicitly. |
| **`AlgSum` for matrix values** | Self-energy Σ = A·I + B·p-slash has no home in scalar `AlgSum` | `AlgSum` is purely commutative scalars. Dirac matrices are non-commutative. | **`DiracExpr`**: separate type for matrix-valued expressions. `DiracExpr = Vector{Tuple{AlgSum, DiracChain}}`. "Everything returns AlgSum" was too aggressive — matrix-valued objects need their own type. |
| **MomSumSlot in traces** | `Tr[γ·(p-k) γ·q]` returned 0 | `gamma_pair` had no MomSumSlot dispatch | **Expand first, trace second.** Don't add complexity to `gamma_pair`. Pre-expand MomentumSum gammas linearly before tracing. |
| **DiracTrick missing** | No gamma contraction within chains (γ^μ...γ_μ) | v2 had trace but not simplification | **Trace ≠ simplification.** These are separate operations. DiracTrick acts on `DiracExpr` (matrix), trace acts on `DiracChain` and returns scalar `AlgSum`. |
| **Spin sum: 1 trace vs 2** | `spin_sum_amplitude_squared` built one big trace instead of Tr₁ × Tr₂ | Two fermion lines need separate traces, not concatenation | **Separate fermion lines produce separate traces.** The product `|M|² = Tr[...] × Tr[...]` is a product of AlgSums, not a single trace. Each line gets its own completeness insertion and trace. The conjugate amplitude reverses gamma order and relabels indices. |
| **Conjugate index relabeling** | Shared `γ^μ` between traces caused self-contraction | Amplitude and conjugate share index names | **Relabel indices in conjugate.** `μ → μ_` (prime) to avoid premature contraction. After tracing, the primed and unprimed indices contract via the photon propagator metric. |
| **DimPoly display** | `show` printed rationals like `-4//1` | Display logic for polynomial terms needs integer detection | Minor. Use `isinteger(c) ? Int(c) : c` in `show`. |
| **AlgFactor union growing** | 6 types: `Pair, Eps, SUNDelta, FundDelta, SUNF, SUND` | Each physics space adds factor types | **Manageable at 6.** Julia handles small unions well. If it grows past ~8, consider an abstract type or two-level union (`LorentzFactor`, `ColourFactor`). |
| **SymDim was too narrow** | `SymDim * SymDim` errored (D² not representable) | Only `aD + b` representable | **DimPoly (polynomial in D)** handles all powers. 2-loop gives D². Cheap to implement, future-proof. |
| **Concrete N for colour** | `N=3` evaluated immediately; symbolic N not available | Design choice for simplicity | **Acceptable for QCD.** For symbolic N (large-N limits, SU(2) comparison), need `NPoly` type analogous to `DimPoly`. Defer until needed. |

### Design decisions that STAYED correct

- **Pair as universal Lorentz bilinear** — never questioned, always right
- **SPContext as explicit context** — ScopedValues adds convenience without losing testability
- **Canonical ordering in constructors** — prevents duplicate terms, enables Dict hashing
- **Copy-on-write for SPContext** — immutable, thread-safe
- **200 LOC file limit** — forced good decomposition

### Key architectural insight

**The coefficient type IS the architecture.** v1's `coeff::Any` + Expr building infected
every function with adaptation code. v2's `DimPoly` eliminates ~150 LOC of glue.
Get the coefficient representation right FIRST, before writing anything else.

Similarly, **the return type IS the API.** "Everything returns AlgSum" eliminated all
normalization layers. The one exception (DiracExpr for matrix values) was forced by physics,
not a design failure.

---

## Type Hierarchy (UPDATED)

```
    PhysicsIndex (abstract)
      ├─ LorentzIndex         (μ, ν, ρ, σ)
      ├─ AdjointIndex         (a, b, c — SU(N) adjoint)
      └─ FundIndex            (i, j, k — SU(N) fundamental)

    PairArg = Union{LorentzIndex, Momentum, MomentumSum}

    Pair{A<:PairArg, B<:PairArg}              # NOT exported
      ├─ MetricTensor  = Pair{LorentzIndex, LorentzIndex}
      ├─ FourVector    = Pair{LorentzIndex, Momentum}
      └─ ScalarProduct = Pair{Momentum, Momentum}

    DiracGamma{S<:DiracSlot}                  # parametric
    Spinor{K<:SpinorKind}                     # parametric

    AlgFactor = Union{Pair, Eps, SUNDelta, FundDelta, SUNF, SUND}
    AlgSum = Dict{FactorKey, Coeff}           # commutative scalars
    DiracExpr = Vector{Tuple{AlgSum, DiracChain}}  # matrix-valued

    Coeff = Union{Rational{Int}, DimPoly}
    DimPoly = polynomial in D                 # handles D, D-4, D², ...

    GaugeGroup: U1, SU{N}                    # compile-time dispatch
    FieldSpecies: Fermion, Boson, Scalar
    Field{S<:FieldSpecies}                    # parametric
    AbstractModel → QEDModel (concrete)
    FeynmanRules (callable struct)
    FeynmanDiagram (callable struct)
    CrossSectionProblem → solve_tree()        # DiffEq pattern
```

---

## The Six-Layer Pipeline

```
Layer 1: Model     →  qed_model()                → QEDModel
Layer 2: Rules     →  feynman_rules(model)        → FeynmanRules (callable)
Layer 3: Diagrams  →  tree_diagrams(model, ...)   → Vector{FeynmanDiagram}
Layer 4: Algebra   →  trace → contract → expand   → AlgSum
Layer 5: Integrals →  (not yet implemented)       → PaVe numerical
Layer 6: Evaluate  →  evaluate_sp → σ_total       → Float64 (nanobarns)
```

Design patterns per layer:
- **Layer 1**: Abstract interface (`AbstractModel`), field species as type param, traits
- **Layer 2**: Callable struct, dispatch on species + gauge group for vertices/propagators
- **Layer 3**: Callable diagram, hard-coded topologies (full generation deferred)
- **Layer 4**: Parametric types, Dict AlgSum, DimPoly coefficients
- **Layer 5**: Native Julia (Li₂, A₀, B₀), no Fortran FFI
- **Layer 6**: DiffEq Problem pattern, Mandelstam kinematics, analytical phase space

---

## Anti-Patterns — DO NOT DO THESE

1. **`x.field isa SomeType`** → parametric type + dispatch
2. **`coeff::Any`** → `Coeff = Union{Rational{Int}, DimPoly}`
3. **Building `Expr` trees** → `DimPoly` arithmetic
4. **Mixed return types** → `AlgSum` for scalars, `DiracExpr` for matrices
5. **`repr()`-based sorting** → structural `Base.isless`
6. **Global mutable state** → `ScopedValues` or explicit context
7. **`Vector{Any}`** → small concrete Unions or parametric containers
8. **One trace for two fermion lines** → separate traces, multiply results
9. **Shared indices between amplitude and conjugate** → relabel conjugate indices
10. **OOP struct hierarchies by default** → research Julia idioms first

---

## How to Add a New Physics Space

Proven pattern (used for colour algebra):

1. Define index types as `<: PhysicsIndex`
2. Define factor types with canonical constructors
3. Add to `AlgFactor` union + implement `isless`, `hash`, `==`
4. Add `_factor_type_tag` for structural ordering
5. `alg(f::AlgFactor)` lifter works automatically
6. Implement trace/simplify returning `AlgSum`
7. Implement contraction as a separate pass

---

## File Map

```
src/v2/
├── FeynfeldX.jl          # Module root (79 LOC)
├── coeff.jl              # DimPoly coefficient algebra (142)
├── types.jl              # PhysicsIndex, DimSlot, LorentzIndex, Momentum (79)
├── colour_types.jl       # AdjointIndex, FundIndex, SUNT, SUNF, SUND, deltas (126)
├── pair.jl               # Parametric Pair{A,B}, not exported (76)
├── expr.jl               # AlgSum (Dict), AlgFactor (6-type union), FactorKey (149)
├── sp_context.jl         # SPContext + ScopedValues (71)
├── contract.jl           # Lorentz index contraction (136)
├── expand_sp.jl          # Scalar product bilinear expansion (82)
├── dirac.jl              # DiracGamma{S}, Spinor{K}, DiracChain (113)
├── dirac_trace.jl        # Dirac trace → AlgSum (65)
├── dirac_expr.jl         # DiracExpr: matrix-valued expressions (104)
├── dirac_trick.jl        # D-dimensional gamma contraction n=0,1,2 (117)
├── spin_sum.jl           # Fermion spin sums: separate traces × multiply (138)
├── colour_trace.jl       # SU(N) trace → AlgSum, concrete N (82)
├── colour_simplify.jl    # Delta contraction, f·f identity (148)
├── model.jl              # AbstractModel, Field{Species}, GaugeGroup types (101)
├── rules.jl              # FeynmanRules callable, dispatch on species (81)
├── diagrams.jl           # FeynmanDiagram callable, hard-coded topologies (75)
├── cross_section.jl      # Mandelstam, Problem pattern, σ_total (92)
├── DESIGN.md             # THIS FILE
└── VERTICAL_PLAN.md      # Full pipeline plan including Stage B (1-loop)

test/v2/
├── test_ee_mumu_x.jl    # Algebra tracer bullet: P&S (5.10) (14 tests)
├── test_coeff.jl         # DimPoly unit tests (29 tests)
├── test_colour.jl        # SU(N) colour algebra (22 tests)
├── test_self_energy.jl   # 1-loop DiracExpr + DiracTrick (25 tests)
└── test_vertical.jl      # FULL PIPELINE: Model→Rules→Diagrams→Algebra→σ (33 tests)

Total: 2,056 source LOC, 630 test LOC, 123 tests, all files < 200 LOC.
```
