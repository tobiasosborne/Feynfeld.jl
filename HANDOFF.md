# HANDOFF — 2026-03-26 (End of Session 3)

## DO NOT DELETE THIS FILE. Read it completely before working.

## TOBIAS'S RULES — FOLLOW TO THE LETTER

1. **SKEPTICISM**: All subagent work, handoffs — verify everything twice.
2. **DEEP BUGS**: Deep, complex, interlocked. Do not underestimate.
3. **NO BANDAIDS**: Best-practices full solutions only.
4. **WORKFLOW**: 3 subagents before any core code change (research source + 2 solutions).
5. **REVIEW**: Rigorous reviewer agent after every core change. No exceptions.
6. **GROUND TRUTH**: Physics is ground truth, not pinned numbers. Tests may be suspect.
7. **TESTING**: Targeted only, or full suite in background.
8. **REPEAT RULES**: Repeat occasionally to maintain focus.
9. **DO NOT UNDERESTIMATE**: This is deeply nontrivial.
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Read-only research/design agents CAN run in parallel.

**NEVER modify TensorGR.jl without explicit permission.** It is an active separate
project with its own workflow and handoff protocol.

---

## Current State — What's Done

**Phase 1 (algebra layer) complete. Tracer bullet (e+e- → mu+mu-) VALIDATED.**

| Phase | Epic | Tasks | Status | Key deliverables |
|-------|------|-------|--------|-----------------|
| 0 | `feynfeld-d4m` | 5/5 | DONE | BMHV dim algebra, Minkowski TensorGR bridge |
| 1a | `feynfeld-43a` | 7/7 | DONE | Pair, Momentum, Contract, ExpandSP, Eps, SPContext |
| 1b | `feynfeld-mpw` | 9/9 | DONE | DiracGamma, DiracTrick, DiracTrace, DiracSimplify, DiracScheme |
| 1c | `feynfeld-sem` | 4/4 | DONE | SUNT, SUNF, SUND, SUNTrace, colour contractions |
| 1d | `feynfeld-9jr` | 6/6 | DONE | PaVe{N}, FAD, Tdec rank 0-2, PaVeReduce (B-functions) |
| 5* | `feynfeld-89y` | 3/6 | PARTIAL | AlgSum expression tree, FermionSpinSum, P&S Eq (5.10) validated |

- **40/80 beads issues closed** (50%)
- **616 tests, ALL PASS**
- **~2,900 source LOC** across 28 Julia files (all under 200 LOC)
- **~2,200 test LOC** including MUnit translations and tracer bullet

### Session 3 additions
- **AlgTerm/AlgSum** (`alg_expr.jl`, `alg_ops.jl`): Expression tree enabling composable algebra. Sum-of-products with arithmetic (+, *, -), contract(AlgSum), expand_scalar_product(AlgSum), evaluate_sp(AlgSum), dirac_trace_alg.
- **FermionSpinSum** (`fermion_spin_sum.jl`): Completeness relations for spin-averaged |M|². Handles massive and massless fermions.
- **e+e- → mu+mu- tracer bullet** (`test/test_ee_mumu.jl`): End-to-end computation validates P&S Eq (5.10). Spin sums → traces → contract → Mandelstam → 8(t²+u²).

### Beads setup note
The Dolt database may not survive across sessions. To restore:
```bash
bd init --force --prefix feynfeld && bd backup restore
```

---

## Architecture — File Map

```
src/
├── Feynfeld.jl                      # Module root, include order matters
├── algebra/
│   ├── dimensions.jl                # Dim4/DimD/DimDm4, dim_contract, dim_trace, to_dim
│   ├── types.jl                     # FeynExpr abstract, LorentzIndex, FourMomentum (legacy), FTerm, Amplitude
│   ├── momentum.jl                  # Momentum, MomentumSum, PairArg union, arithmetic
│   ├── pair.jl                      # Pair (universal bilinear), pair(), SP/FV/MT + D/E variants
│   ├── sp_context.jl                # SPContext (copy-on-write scalar product storage)
│   ├── expand_sp.jl                 # expand_scalar_product() bilinear expansion
│   ├── eps.jl                       # Eps (Levi-Civita), levi_civita(), LC/LCD, eps_contract()
│   ├── contract.jl                  # contract() worklist Lorentz index contraction (full BMHV)
│   ├── minkowski.jl                 # minkowski_registry() TensorGR bridge
│   ├── dirac_types.jl               # DiracGamma (LISlot/MomSlot/SpecialSlot), Spinor, DiracChain
│   ├── dirac_chain.jl               # dot() constructor, dot_simplify(), DiracElement
│   ├── dirac_trick.jl               # dirac_trick() — core γ-matrix simplification rules
│   ├── dirac_trace.jl               # dirac_trace() — recursive trace formula + γ5 chiral base
│   ├── dirac_order.jl               # dirac_order() — normal ordering via anticommutation
│   ├── dirac_equation.jl            # dirac_equation() — p-slash u(p) = m u(p) at boundaries
│   ├── dirac_simplify.jl            # dirac_simplify() — master orchestrator
│   ├── dirac_scheme.jl              # DiracScheme enum (NDR/BMHV/Larin), with_scheme()
│   ├── colour_types.jl              # SUNIndex/SUNFIndex, SUNT, SUNTF, SUNF, SUND, deltas, ColourChain
│   ├── colour_trace.jl              # SUNTrace, sun_trace() recursive
│   ├── colour_simplify.jl           # delta_trace, contract_ff/dd/fd, structure constant identities
│   ├── alg_expr.jl                  # AlgTerm/AlgSum types, arithmetic (+, *, -)
│   ├── alg_ops.jl                   # contract/expand_sp/evaluate_sp(AlgSum), dirac_trace_alg
│   └── fermion_spin_sum.jl          # FermionSpinSum: completeness relations for |M|²
└── integrals/
    ├── pave.jl                      # PaVe{N} parametric, A0/A00/B0/B1/B00/B11/C0/D0
    ├── feynamp_denominator.jl       # FeynAmpDenominator, 3 propagator types, FAD()
    ├── tdec.jl                      # tdec() tensor decomposition rank 0/1/2
    └── pave_reduce.jl               # pave_reduce() B-function coefficient reduction
```

### Key design decisions (read before changing anything)
- **Pair is the universal Lorentz bilinear** (not exported — conflicts with Base.Pair). Users use `SP`/`FV`/`MT` convenience or `Feynfeld.Pair`. This mirrors FeynCalc's `Pair[x,y]`.
- **BMHV projection fires at Pair construction time** via `dim_contract()`. The `pair()` factory returns `0` on vanishing (4 ∩ (D-4) = 0).
- **SPContext is explicit, not global state.** Passed as `ctx` kwarg to `contract()` and `expand_scalar_product()`.
- **DiracGamma uses tagged DiracSlot union** (LISlot, MomSlot, SpecialSlot) — not separate types.
- **PaVe{N}** is parametric on number of propagators for dispatch in reduction.
- **LorentzIndex default dim is Dim4()** (matches FeynCalc convention).
- **FeynRules port MUST use second-quantization operator algebra** (Tobias rejected functional differentiation).

---

## Known Limitations — What's Incomplete in Phase 1

These are areas where the Phase 1 implementation is deliberately minimal. Future phases may need to extend them.

| # | Area | Gap | Impact |
|---|------|-----|--------|
| 1 | **Symbolic coefficients** | `_mul_coeff` produces raw `Expr` trees (`:($a * $b)`) that don't simplify. `:(D * 4)` stays as-is, not `:(4D)`. | Will bite in Dirac trace results and cross-sections. Need a minimal symbolic simplifier before Phase 5. |
| 2 | **Expression tree** | ~~FIXED in Session 3~~. AlgTerm/AlgSum provides sum-of-products with full arithmetic and contract/expand_sp/evaluate_sp integration. | Used for tracer bullet. FTerm/Amplitude still unused (reserved for Feynman amplitudes). |
| 3 | **DiracTrick sandwich n≥3** | `g^μ g^a g^b g^c g_μ` (3+ gammas between) not implemented. General formula exists (Mertig/Boehm/Denner 1991 Eq 2.9). | Needed for any 1-loop QED/QCD calculation with ≥3 internal propagators. |
| 4 | **DiracTrace chiral n>4** | `Tr[γ5 g^a g^b g^c g^d g^e g^f]` (6+ gammas with γ5) uses recursive reduction — not implemented. | Needed for axial-vector processes. |
| 5 | **BMHV/Larin scheme** | Enum registered, `with_scheme()` works, but DiracTrick doesn't dispatch on scheme yet. All γ5 algebra uses NDR (anticommuting γ5). | BMHV needed for correct D-dim γ5 in 2-loop calculations. |
| 6 | **SUNSimplify Cvitanovic** | Fierz identity `T^a_{ij} T^a_{kl} = (1/2)δ_{il}δ_{kj} - (1/(2N))δ_{ij}δ_{kl}` not implemented. Only delta_trace and f/d contractions. | Needed for any multi-loop colour factor calculation. |
| 7 | **Tdec rank ≥3** | Only rank 0/1/2 tensor decomposition. Rank 3+ needs Passarino-Veltman projection (solving linear system) or TIDL lookup. | Rank 3 appears in box diagrams. |
| 8 | **PaVeReduce C/D** | Only B1/B00/B11 → A0+B0. Denner recursion for C/D coefficients not implemented. | Needed for 3- and 4-point processes. |
| 9 | **ToPaVe** | No FAD → PaVe conversion. Need Denner kinematic invariant extraction from propagator momenta. | Required to connect diagram output to scalar integrals. |
| 10 | **PaVeAutoOrder** | PaVe sorts indices but doesn't canonicalize invariant/mass ordering under permutation symmetries. | May cause missed simplifications in expressions with multiple PaVe symbols. |
| 11 | **GenPaVe** | Not implemented. Needed when physical momentum routing differs from Denner convention. | Required for ToPaVe with non-standard routing. |
| 12 | **No `Base.show` methods** | All types print as raw struct output. No `p^μ` or `g^{μν}` rendering. | Poor ergonomics for interactive use. |

---

## What To Do Next — Recommended Order

**Tracer bullet (Option A) is DONE.** e+e- → mu+mu- tree-level validated.

### Recommended next: extend the pipeline

1. **Phase 2** (`feynfeld-62k`): LoopTools numerical integration
   - Start with Li2 (dilogarithm) — atomic building block
   - Then A0, B0 numerical evaluation
   - Reference: `refs/LoopTools/` has the Fortran source
2. **Phase 3** (`feynfeld-ntj`): Model + Rules (FeynRules port)
   - Must use second-quantization operator algebra (Tobias's requirement)
3. **Phase 4** (`feynfeld-c0n`): Diagrams (FeynArts port)
4. **Phase 5 remaining** (`feynfeld-89y`): PolarizationSum, ColourME, phase-space integration
5. **Phase 6** (`feynfeld-42d`): ULDM application

### Algebra gaps to fix (for 1-loop calculations)
- **Limitation #1** (symbolic simplification): `_mul_coeff` produces raw Expr trees. Need simplifier.
- **Limitation #3** (DiracTrick sandwich n≥3): Needed for 1-loop QED/QCD.
- **Limitation #6** (Fierz identity): Needed for multi-loop colour factors.
- **Limitation #7** (Tdec rank ≥3): Needed for box diagrams.

### Cross-cutting tasks (do anytime)
- `feynfeld-lc0`: FCI/FCE convenience API (macros/constructors) — nice-to-have
- `feynfeld-zwc`: Expression simplify() pipeline — would fix limitation #1
- `feynfeld-2wx`: Feynman rules (propagators, vertices) — bridges Model→Diagrams
- Limitation #12: Base.show methods for interactive use

---

## Reference Codebases

All in `refs/` (gitignored):
- `refs/FeynCalc/` — Primary porting oracle. 186k LOC Mathematica. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/FormCalc/` — FormCalc reference (NOT directly ported).
- `refs/LoopTools/` — Loop integral numerics: dual FF+Denner Fortran implementation.

Architecture reports from Session 1 in `refs/reports/`:
- `feyncalc_architecture.md` — Full FeynCalc deep-dive
- `feynarts_architecture.md` — FeynArts deep-dive
- `feynrules_architecture.md` — FeynRules deep-dive
- `formcalc_looptools_architecture.md` — FormCalc+LoopTools deep-dive
- `tensorgr_patterns.md` — TensorGR integration patterns

---

## Key Dependencies

- **TensorGR.jl** (`../TensorGR.jl`): Shared index contraction/canonicalization engine.
  Feynfeld's `minkowski_registry()` uses it. DO NOT modify without permission.
- **Julia 1.12+**: Package uses modern Julia features.
- **No external dependencies** beyond TensorGR and stdlib (LinearAlgebra, Test).

---

## Quick Commands

```bash
# Build & test
julia --project=. -e 'using Pkg; Pkg.test()'     # 506 tests

# Beads issue tracking
bd init --force --prefix feynfeld && bd backup restore  # restore after session break
bd stats                                            # 36 closed, 41 open
bd ready                                            # available work
bd show feynfeld-62k                                # Phase 2 epic
bd show feynfeld-89y                                # Phase 5 epic
bd list --status=open --limit 0                     # all open issues

# Git
git log --oneline -10                               # recent commits
```
