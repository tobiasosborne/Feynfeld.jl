# CLAUDE.md — Feynfeld.jl

## What is this?
Feynfeld.jl is a Julia package for symbolic and numerical quantum field theory computation.
It replaces the Mathematica ecosystem (FeynRules, FeynArts, FeynCalc, FormCalc, LoopTools)
with a single end-to-end pipeline: Lagrangian → cross-section in one `using Feynfeld`.

## Architecture: six layers
```
Model → Rules → Diagrams → Algebra (core) → Integrals → Evaluate
```
The Algebra layer is the type system. Build order follows the PRD (Phase 0–6).
See `Feynfeld_PRD.md` for full details.

## Key dependencies
- **TensorGR.jl** (`../TensorGR.jl`): Shared index contraction/canonicalization engine.
  Feynfeld's Lorentz algebra specialises TensorGR's abstract tensor layer to η^{μν}.

## Reference codebases (ground truth)
All in `refs/` (gitignored):
- `refs/FeynCalc/` — Primary porting oracle. MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/FormCalc/` — FormCalc reference (not directly ported).
- `refs/LoopTools/` — Loop integral numerics reference.

## Build & test
```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## TOBIAS'S RULES — FOLLOW TO THE LETTER
1. **SKEPTICISM**: Verify subagent work twice.
2. **DEEP BUGS**: Do not underestimate complexity.
3. **NO BANDAIDS**: Best-practices full solutions only.
4. **WORKFLOW**: 3 subagents before any core code change (research source + 2 solutions).
5. **REVIEW**: Rigorous reviewer agent after every core change.
6. **GROUND TRUTH**: Physics is ground truth, not pinned numbers. Local copies of published
   papers or reference implementations are the ONLY truth.
7. **TESTING**: Targeted only, or full suite in background.
8. **DO NOT UNDERESTIMATE**: This is deeply nontrivial.
9. **NO PARALLEL AGENTS**: Julia precompilation cache conflicts. Run agents sequentially.
10. **LOC LIMIT**: No source file exceeds ~200 lines. Split aggressively.

## MUnit translation protocol
For each FeynCalc MUnit test file in `refs/FeynCalc/Tests/`:
1. Read the `.mt` file.
2. Translate each test to a Julia `@test`, preserving math exactly.
3. Document the source MUnit file and test ID in a comment.
4. If MUnit test and textbook disagree, the textbook wins (Rule 6).
