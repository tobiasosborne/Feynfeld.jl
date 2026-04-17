# Feynfeld.jl v2 Algebra Layer (Layer 4) Stocktake — 2026-04-17

**Overall:** 20 files, ~2,040 LOC. Core algebra is **excellent** (validated by Session 8 review). Type system is Julian, dispatch-based, type-stable except for known issues in factory functions.

---

## File-by-File Survey

| File | LOC | Purpose | Public API | Status |
|------|-----|---------|-----------|--------|
| **coeff.jl** | 142 | DimPoly polynomial in D for dim-reg coefficients | DimPoly, DIM, DIM_MINUS_4, mul_coeff, add_coeff | Solid |
| **types.jl** | 101 | Core: DimSlot, LorentzIndex, Momentum, MomentumSum | Dim4, DimD, DimDm4, LorentzIndex, Momentum, momentum_sum factory | Good (factory returns Union) |
| **pair.jl** | 75 | Parametric Pair{A,B} Lorentz bilinear + user API | pair(), SP(), FV(), MT(), SPD(), FVD(), MTD() | Excellent |
| **expr.jl** | 150 | AlgSum Dict-based expression (FactorKey → Coeff) | AlgSum, alg() constructor, Eps, AlgFactor union | Excellent |
| **contract.jl** | 212 | Einstein index contraction via dispatch on Pair types | contract(), substitute_index() | Excellent |
| **expand_sp.jl** | 95 | Bilinear expansion of MomentumSum in Pair/Eps | expand_scalar_product() | Good (Tuple{Any,...} issue) |
| **dirac.jl** | 120 | DiracGamma{S} parametric type, Spinor{K}, DiracChain | GA, GS, GA5, GA6, GA7, u, v, ubar, vbar, dot() | Excellent |
| **dirac_expr.jl** | 107 | DiracExpr matrix-valued Dirac algebra (scalar*chain pairs) | DiracExpr, arithmetic, simplify() | Good |
| **dirac_trace.jl** | 158 | Dirac trace with γ5 expansion → AlgSum | dirac_trace(), returns uniform AlgSum | Excellent |
| **dirac_trick.jl** | 186 | γ^μ Γ^(l) γ_μ sandwich contraction (Eq. 2.9) | dirac_trick() | Excellent |
| **eps_contract.jl** | 86 | ε·ε → determinant of metric tensors (S₄ perms precomputed) | eps_contract() | Good |
| **eps_evaluate.jl** | 65 | Numerical eval of fully-contracted ε(p,q,r,s) via Gram det | _evaluate_eps(), Gram-determinant sign convention | Good |
| **colour_types.jl** | 122 | SU(N) colour algebra: AdjointIndex, FundIndex, SUNDelta, SUNF, SUND | Parametric types, _fresh_adj() gensym | Excellent |
| **colour_simplify.jl** | 183 | Delta contraction, structure constant identities (f·f, d·d, f·d) | contract_colour() | Excellent |
| **colour_trace.jl** | 84 | Colour trace with real/imag splitting for i² handling | colour_trace(), Tr(T^a T^b T^c) = (1/4)(d+if) | Excellent |
| **spin_sum.jl** | 171 | Fermion completeness + spin-summed |M|² (two separate traces) | fermion_spin_sum(), spin_sum_amplitude_squared() | Good (imports from polarization_sum) |
| **polarization_sum.jl** | 57 | Gauge boson polarization: Feynman + axial + massive | polarization_sum(), polarization_sum_massive() | Good |
| **interference.jl** | 106 | Cross-line traces for multi-channel interference (Bhabha s×t) | spin_sum_interference() | Good |
| **sp_context.jl** | 71 | Scalar product context via ScopedValues (Julia 1.11+) | SPContext, with_sp(), set_sp(), sp_context() | Good |
| **sp_lookup.jl** | 32 | Float64 SP lookup + K·p dot products for TID | _Kdot(), _splookup() | Minimal |

---

## Cross-cutting Observations

### Type Hierarchy & Dispatch
- **PairArg union** (LorentzIndex, Momentum, MomentumSum) enables dispatch on Pair type parameters:
  - `Pair{LI, LI}` → metric tensor
  - `Pair{LI, M}` → four-vector
  - `Pair{M, M}` → scalar product
  - Clean parametric design avoids isa cascades.

- **AlgFactor union** (6 types: Pair, Eps, SUNDelta, FundDelta, SUNF, SUND) stored in FactorKey → enables structural ordering via dispatch table `_factor_type_tag()`.

- **Dirac slot dispatch** (DiracGamma{S} where S <: DiracSlot) replaces symbol checking: LISlot, MomSlot, Gamma5Slot, ProjPSlot, etc. **Type-stable pattern.**

- **Spinor{K} parametric types** (UKind, VKind, UBarKind, VBarKind) dispatch on completeness relations without symbol matching.

### Known Type Instabilities (CLAUDE.md, Session 8)
1. **momentum_sum() factory** returns `Union{Nothing, Momentum, MomentumSum}` (line 62–71, types.jl)
2. **gamma_pair() fallback** returns `Union{Int, Pair, AlgSum}` (dirac.jl:120, maps 0 → AlgSum())
3. **pair() factory** returns `Union{Int, Pair}` when BMHV vanishing occurs (returns 0 vs Pair)
4. **expand_sp.jl lines 41, 54, etc.** use `Tuple{Coeff, Union{AlgFactor, Nothing}}` (could be typed as Coeff ⊔ AlgFactor)
5. **Global mutable:** `_fresh_adj()` uses gensym (good), **but** colour_simplify.jl references (no global counter detected — clean)

### Coefficient System
- **Coeff = Union{Rational{Int}, DimPoly}** (coeff.jl:121) — maximally efficient for Julia union-splitting.
- **normalise_coeff()** collapses DimPoly constants back to Rational{Int} — memory-efficient.
- **DimPoly invariant:** no trailing zeros; empty = zero. Constructor enforces.
- **mul_coeff() / add_coeff()** dispatch cleanly — no Any anywhere.

### Contraction & Expansion Strategy
- **contract()** (contract.jl:4–12): two-pass architecture:
  1. eps_contract() pre-pass (determinant formula, S₄ permutation table)
  2. Per-index worklist loop (repeated scans for contractable pairs)
  
- **expand_scalar_product()** (expand_sp.jl:4–31): recursive expansion of MomentumSum in Pair/Eps slots. Returns AlgSum always.

- **Index substitution** (contract.jl:181–205) for polarization sum relabeling — clean dispatch on Pair type.

### Dirac Algebra
- **dirac_trace()** handles MomentumSum expansion first (linearity), then projects/γ5 via **recursive pre-processing** (lines 19–46, dirac_trace.jl):
  - MomentumSum expands by coefficient.
  - Projectors expand to (I ± γ5)/2.
  - γ5 anticommutes to end, sign tracked.
  - **Returns AlgSum always** — no mixed return types.

- **dirac_trick()** (dirac_trick.jl:16–186): sandwich formula γ^μ Γ^(l) γ_μ dispatches on inner-gamma count (n=0,1,2,3,4,≥5). Refs Mertig1991 Eq. 2.9. **Verified against FeynCalc tests.**

- **DiracExpr simplify()** (dirac_expr.jl:93–107): groups terms by gamma structure (using gammas() vector as key). Combines identical chains. **Not lazy — eager simplification after each operation.**

### Colour Algebra
- **Parametric colour indices** (AdjointIndex, FundIndex) mimic LorentzIndex pattern.
- **SUNF/SUND canonical sorting with sign tracking** — antisymmetric (SUNF) uses parity of permutation; symmetric (SUND) sorts unconditionally.
- **colour_trace()** tracks real + imaginary parts separately (lines 26–68, colour_trace.jl) to handle i² = -1 correctly for f-term products. **Returns real part only** (sufficient for |M|²).
- **structure constant identities** (colour_simplify.jl:132–163):
  - f·f: 3-index match → N(N²-1), 2-index → N δ^ab
  - d·d: symmetry-respecting contraction
  - f·d: vanishes by index antisymmetry

### Epsilon Tensor (Levi-Civita)
- **eps_contract()** pre-computes S₄ permutation table (24 perms with parity). **O(1) per Eps pair.**
- **Determinant formula** (Mertig1991 Eq. 2.21): ε₁·ε₂ = -det[pair(a_i, b_j)]_{4×4}. Uses Leibniz summation over permutations.
- **eps_evaluate()** computes |ε| from Gram determinant (`√(-det G)`), sign via reference order (:p1,:p2,:k1,:k2). **Convention: ε > 0 when θ ∈ (0,π) in CM frame.**

### Spin Sums & Polarization
- **spin_sum_amplitude_squared()** (spin_sum.jl:64–73): **TWO separate traces multiplied**, not one big trace. Key insight for amplitude-squared computation.
- **_single_line_trace()** (spin_sum.jl:77–96): DiracExpr case loops over all (i,j) term pairs, applies conjugate relabeling to avoid index collision.
- **_conjugate_gammas()** (spin_sum.jl:127–140): reverses gamma order AND relabels indices (μ → μ_). **Projector conjugation:** GA6 ↔ GA7.
- **polarization_sum()** three forms:
  1. Feynman (massless): -g^{μν}
  2. Axial (with reference vector n): -g^{μν} + (k^μ n^ν + n^μ k^ν)/(k·n)
  3. Massive: -g^{μν} + k^μ k^ν / M²

### Scalar Product Context
- **SPContext** uses Julia 1.11+ ScopedValues for implicit threading. **Not global mutable state — functional.**
- **sp_lookup.jl** minimal helpers for TID evaluation (Gram matrices, K·p dot products).

---

## Tech Debt & Known Issues

### Type Instability
- **momentum_sum(), gamma_pair(), pair() factories** return unions. See CLAUDE.md Session 8 for remediation priority.
- **expand_sp.jl result tuples** use `Tuple{Coeff, Union{AlgFactor, Nothing}}` — consider type alias for clarity.

### Missing Functionality
- **gamma5 not fully implemented** — dirac_trace.jl has full γ5 algebra, but **no algebraic simplification of γ5 in dirac_trick.jl** yet (Spiral 8 task).
- **Eps slot with MomentumSum in dirac_trace** expanded (line 38–30), but **not in dirac_trick** (would need parametric dispatch on MomSumSlot in _trick_gammas).

### Coupling & Layering
- **spin_sum.jl imports from polarization_sum.jl** (implicit via _lookup_sp). No circular dependency but tightly coupled.
- **interference.jl** copies _completeness, _conjugate_gammas from spin_sum.jl (code duplication, not modularized).
- **colour_simplify.jl line 51:** `_fresh_adj()` uses gensym — clean, but not thread-safe for concurrent simplifications (ScopedValue pattern would be safer).

### Documentation & Testing Gaps
- No TODO/FIXME in algebra layer (good). But CLAUDE.md MUST FIX list (type instabilities) is not yet targeted.
- **Cross-references to FeynCalc tests present** (e.g., dirac_trick.jl:121, eps_contract.jl:7), but **no explicit test file references in v2/algebra/** — tests live in `test/v2/` separately.

### Experimental Code
- **d_tensor.jl, degree_partition.jl, d0_collier.jl** (loop integral helpers) exist but **not in this stocktake scope** (Layer 5: Integrals).
- **ew_model.jl, ew_cross_section.jl, ew_parameters.jl** (EW algebra) — **not surveyed** (Model/Evaluate layers).

---

## Summary Assessment

| Aspect | Rating | Notes |
|--------|--------|-------|
| **Type Design** | Excellent | Parametric, dispatch-based, Julian. Avoids isa cascades. |
| **Coefficient System** | Excellent | DimPoly with proper normalization. Union{Rational, DimPoly} efficient. |
| **Lorentz Algebra** | Excellent | Pair parametricity, contract/expand dispatch clean. |
| **Dirac Algebra** | Very Good | Traces uniform AlgSum. γ5 recursion correct. Trick formula verified. |
| **Colour Algebra** | Excellent | Canonical forms, sign tracking, trace with i² handling. |
| **Epsilon Tensor** | Very Good | S₄ precomputed, determinant formula clean. Sign convention clear. |
| **Spin Sums** | Very Good | Two-trace architecture correct, conjugate relabeling works. |
| **Factory Stability** | Poor | momentum_sum, gamma_pair, pair return unions. **MUST FIX (Session 8).** |
| **Code Coupling** | Good | Minimal cross-imports. interference.jl duplication only issue. |
| **Completeness** | Very Good | No TODO in algebra. γ5 simplification deferred (Spiral 8). |

**Next Actions:** Fix type instabilities in factories (pair.jl, dirac.jl, types.jl) → enable safe dispatch in downstream layers.

