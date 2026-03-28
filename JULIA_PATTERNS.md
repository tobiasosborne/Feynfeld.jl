# Julia Idiom Cheatsheet for Feynfeld.jl

**READ THIS BEFORE WRITING ANY CODE.** This is replicated in `Feynfeld_PRD.md` §6
and `CLAUDE.md`. All three copies must stay in sync.

## DO

```julia
# Parametric types for dispatch (not tagged unions)
struct DiracGamma{S<:DiracSlot}
    slot::S
end

# Multiple dispatch (not if/elseif isa cascades)
contract(a::MetricTensor, b::FourVector, idx) = ...
contract(a::FourVector, b::FourVector, idx) = ...

# Dict-based expression storage (O(1) like-term collection)
struct AlgSum
    terms::Dict{FactorKey, Coeff}
end

# Concrete small Unions (Julia optimises up to ~4 types)
const Coeff = Union{Rational{Int}, DimPoly}

# Named constructors (plain functions returning types)
A0(m2::Real) = PaVe{1}(Int[], Float64[], [Float64(m2)])

# ScopedValues for implicit context
const CURRENT_SP = ScopedValue(SPContext())

# Use existing packages
using QuadGK: quadgk
using PolyLog: li2

# evaluate via dispatch, not type-in-function-name
evaluate(pv::PaVe{1}; mu2=1.0) = ...
evaluate(pv::PaVe{2}; mu2=1.0) = ...
```

## DO NOT

```julia
# ✗ Tagged union with isa checks
if x isa ScalarProduct ... elseif x isa MetricTensor ... end

# ✗ Any or Expr as coefficient type
struct BadTerm; coeff::Any; end

# ✗ Building Expr trees
coeff = :($(a) * $(b) + $(c))

# ✗ Type-in-function-name (Java pattern)
evaluate_pave(x)     # ✗
evaluate(x::PaVe)    # ✓

# ✗ Wrapper structs when a function will do
struct SchwingerProblem; alpha::Float64; end   # ✗
schwinger_correction(; alpha=1/137.036) = ...  # ✓

# ✗ Hand-rolled numerical methods
const MY_GAUSS_NODES = (...)  # ✗
quadgk(f, 0, 1)              # ✓

# ✗ Reimplementing standard functions
function my_dilog(z) ... end  # ✗
using PolyLog: li2            # ✓

# ✗ OOP struct hierarchies by default
abstract type AbstractIntegral end  # ✗ (unless dispatch genuinely needed)

# ✗ Forcing types into unions they don't belong in
const AlgFactor = Union{Pair, Eps, ..., PaVe}  # ✗ PaVe is scalar, not a tensor factor

# ✗ Global mutable state
CURRENT_MODEL = nothing  # ✗
const CURRENT_MODEL = ScopedValue(nothing)  # ✓
```

## Julia ecosystem patterns to follow

| Pattern | Example package | What to learn |
|---------|----------------|---------------|
| Dict-based expression storage | SymbolicUtils.jl | `AlgSum` pattern |
| Problem/Solve | DiffEq.jl | `CrossSectionProblem` → `solve_tree` |
| Abstract interface + traits | MultivariatePolynomials.jl | `AbstractModel` pattern |
| Callable structs | Flux.jl | `FeynmanRules(fields)` |
| ScopedValues | Julia 1.11+ stdlib | `SPContext` pattern |
