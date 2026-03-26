# HANDOFF — 2026-03-26 (Session 1)

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
10. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Run agents sequentially only.

**NEVER modify TensorGR.jl without explicit permission.** It is an active separate
project with its own workflow and handoff protocol.

---

## Current State

- **Phase 0 in progress** (TensorGR bridge + BMHV dimension algebra)
- **Full test suite: 63 tests, ALL PASS**
- Beads issue tracker initialized: `bd list` for all issues, `bd ready` for available work
- 7 epics, 67 open tasks across Phases 0–6 (see `bd list --limit 0`)

---

## What Was Done This Session

### 1. Project setup on new machine

- Cloned 5 reference codebases into `refs/` (gitignored):
  - `refs/FeynCalc/` from github.com/FeynCalc/feyncalc
  - `refs/FeynArts/` from github.com/FeynCalc/feynarts-mirror
  - `refs/FeynRules/` from github.com/FeynRules/feynrules
  - `refs/FormCalc/` from github.com/HEPcodes/FormCalc
  - `refs/LoopTools/` from github.com/HEPcodes/LoopTools
- Dev-linked TensorGR.jl as local dependency (`Pkg.develop(path="../TensorGR.jl")`)
- Initialized beads issue tracker (`bd init --prefix feynfeld`)

### 2. Deep architecture research (5 Opus agents, reports in refs/reports/)

Each agent produced a detailed markdown report:

| Report | File | Key findings |
|--------|------|-------------|
| **FeynCalc** | `refs/reports/feyncalc_architecture.md` | `Pair[x,y]` is the universal building block. 97k LOC, 18k MUnit tests. FCI/FCE dual representation. PairContract UpValues trick for contraction. BMHV dim algebra baked into constructors. |
| **FeynArts** | `refs/reports/feynarts_architecture.md` | 3-level Generic→Classes→Particles model hierarchy. Recursive topology generation. Core physics ~2800 lines. |
| **FeynRules** | `refs/reports/feynrules_architecture.md` | Second-quantization operator algebra for vertex extraction (NOT functional differentiation — Tobias rejected that approach). Core path ~7600 lines. |
| **FormCalc+LoopTools** | `refs/reports/formcalc_looptools_architecture.md` | FormCalc NOT directly ported. LoopTools IS porting target — dual FF+Denner implementation, Li2 is atomic building block. |
| **TensorGR.jl** | `refs/reports/tensorgr_patterns.md` | Flat 5-type AST, TIndex with vbundle, Butler-Portugal canonicalization. Key: use TIndex directly, register Minkowski metric, delegate contraction. |

### 3. Granular implementation plan (77 beads issues)

7 epics across 6 phases:
- **Phase 0** (`feynfeld-d4m`): TensorGR bridge — 5 tasks
- **Phase 1a** (`feynfeld-43a`): Lorentz algebra (Pair/Contract) — 7 tasks
- **Phase 1b** (`feynfeld-mpw`): Dirac algebra — 9 tasks
- **Phase 1c** (`feynfeld-sem`): SU(N) colour algebra — 4 tasks
- **Phase 1d** (`feynfeld-9jr`): PaVe/tensor decomposition — 6 tasks
- **Phase 2** (`feynfeld-62k`): LoopTools numerical integrals — 9 tasks
- **Phase 3** (`feynfeld-ntj`): Model + Rules (FeynRules port) — 9 tasks
- **Phase 4** (`feynfeld-c0n`): Diagrams (FeynArts port) — 6 tasks
- **Phase 5** (`feynfeld-89y`): Evaluate (amplitude squaring) — 6 tasks
- **Phase 6** (`feynfeld-42d`): ULDM application — 3 tasks
- Cross-cutting: 3 tasks

19 tasks labelled `so` (spike/research) as placeholders for research rounds.

### 4. Phase 0 implementation (partial)

#### TensorGR.jl changes (committed as d84aa00, pushed to master)
- Widened `ManifoldProperties.dim` and `VBundleProperties.dim` from `Int` to `Union{Int,Symbol}`
- Guarded `contract_metrics` metric trace, `define_metric!` epsilon/signature, `_fs_ddi_order` DDI capping
- All 375,404 TensorGR tests pass, Tier 1 benchmarks pass
- Full details in TensorGR.jl HANDOFF.md

#### Feynfeld.jl new files
- `src/algebra/dimensions.jl`: BMHV dimension algebra — `Dim4`, `DimD`, `DimDm4` types with `dim_contract` dispatch table (9 projection rules), `dim_trace`, `to_dim`
- `src/algebra/minkowski.jl`: `minkowski_registry()` — creates TensorGR registry with Minkowski manifold, flat metric η, supports `dim=4` and `dim=:D`
- `test/algebra/test_dimensions.jl`: 23 tests for BMHV algebra
- `test/algebra/test_minkowski.jl`: 12 tests for Minkowski registry (metric contraction, dim trace)

#### Architectural decision: LorentzIndex stays Feynfeld-specific
FeynCalc's BMHV scheme puts dimension info at the **index level** (individual indices can be 4D, D-dim, or (D-4)-dim). TensorGR puts dimension at the **manifold level**. These don't map cleanly. Decision: keep Feynfeld's `LorentzIndex` struct (already matches FeynCalc semantics), bridge to TensorGR via `minkowski_registry()` when contraction/canonicalization is needed.

---

## Key Decisions / Lessons

- **FeynRules port MUST follow reference implementation**: Second-quantization operator algebra, NOT functional differentiation. Tobias explicitly rejected the alternative.
- **NEVER modify TensorGR.jl without explicit permission**: It's a separate active project. Always check repo status, read HANDOFF.md, ask first.
- **LorentzIndex is Feynfeld-specific, not a TIndex alias**: BMHV dimension tagging at index level doesn't map to TensorGR's manifold-level dimensions.
- **FCI/FCE unnecessary in Julia**: Constructors ARE the internal form. Provide convenience macros instead.
- **Research agents for Mathematica codebases can run in parallel**: Rule 10 (no parallel agents) is about Julia precompilation conflicts. Read-only Mathematica research has no such issue.

---

## ⚠ Core Changes To Monitor

**Phase 0 code** (dimensions.jl, minkowski.jl):
- Location: `src/algebra/dimensions.jl`, `src/algebra/minkowski.jl`
- Risk: Low — foundation types, no existing code depends on them yet
- 63 tests pass

**TensorGR.jl dependency** (commit d84aa00):
- Feynfeld depends on the symbolic-dim support added this session
- If TensorGR.jl reverts the `Union{Int,Symbol}` change, Feynfeld's `minkowski_registry(dim=:D)` will break
- TensorGR.jl TODO: add dedicated symbolic-dim tests, update docstrings

---

## TODO Next Session

1. **Close Phase 0 beads issues** (`feynfeld-d4m.1` through `feynfeld-d4m.5`)
2. **Begin Phase 1a**: Design Pair type system (`feynfeld-43a.1` — spike/research)
3. **Implement Momentum type** (`feynfeld-43a.2`)
4. **Implement ScalarProduct** (`feynfeld-43a.3`)
5. **Implement ExpandScalarProduct** (`feynfeld-43a.4`)
6. **Implement Contract** (`feynfeld-43a.5`) — the big one

## Quick Commands

```bash
bd list --limit 0              # all issues
bd ready                       # available work
bd show feynfeld-d4m           # Phase 0 epic
julia --project=. -e 'using Pkg; Pkg.test()'  # full test suite (63 tests)
```
