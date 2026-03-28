# CLAUDE.md — Feynfeld.jl

## What is this?
Feynfeld.jl is a Julia-native, agent-facing, full-stack physics computation suite.
Lagrangian → cross-section in one `using Feynfeld`. Replaces the Mathematica ecosystem
(FeynRules, FeynArts, FeynCalc, FormCalc, LoopTools) with a single Julia package.

See `Feynfeld_PRD.md` for full vision, architecture, and spiral development plan.
See `src/v2/DESIGN.md` for v2 design choices and anti-patterns.

## Architecture: six layers
```
Model → Rules → Diagrams → Algebra (core) → Integrals → Evaluate
```
The Algebra layer is the type system. Development follows the spiral methodology
(PRD §3): each process drives both vertical (pipeline) and horizontal (MUnit test
coverage) progress.

## Active code
- **v2** (`src/v2/`, branch `experimental/rebuild-v2`): Active development. 187 tests.
- **v1** (`src/algebra/`, `src/integrals/`): FROZEN. Will be deleted. Do not extend.
  Do not import patterns from. v1 used Mathematica-shaped anti-patterns.

## Reference codebases (ground truth)
All in `refs/` (gitignored):
- `refs/FeynCalc/` — Primary porting oracle. 15,222 MUnit tests in `Tests/`.
- `refs/FeynArts/` — Diagram generation reference.
- `refs/FeynRules/` — Model/Lagrangian reference.
- `refs/LoopTools/` — Loop integral numerics (Fortran source).
- `refs/papers/` — Local copies of published papers (Denner 1993, P&S, etc.).

## Ground truth acquisition

Every physics formula must cite a locally stored source (Rule 2). If a paper is
not in `refs/papers/`, acquire it BEFORE writing any code that depends on it.

**Retrieval priority:**
1. **arXiv** — free, use `WebFetch` directly
2. **TIB VPN** — Tobias has institutional access at LUH. Ask him to connect if
   paywalled papers are needed (APS, Springer, Elsevier).
3. **Playwright-cli + headed Chrome** — for paywalled PDFs behind Cloudflare or
   institutional auth. See `../FQHE/scripts/fetch_via_browser.sh` for the pattern:
   ```bash
   # Requires: TIB VPN active, browser session authenticated
   playwright-cli run-code "async page => {
     const resp = await page.request.get('URL', { timeout: 30000 });
     const body = await resp.body();
     return body.toString('base64');
   }"
   ```
   This uses the existing browser session cookies to bypass paywalls.

**Required papers (minimum set):**
- Peskin & Schroeder (1995) — `refs/papers/` (present, .djvu format)
- Denner, Fortschr. Phys. 41 (1993) 307 — one-loop techniques, PV reduction
- 't Hooft & Veltman, Nucl. Phys. B153 (1979) 365 — scalar integrals
- Passarino & Veltman, Nucl. Phys. B160 (1979) 151 — PV decomposition
- Shtabovenko et al., "FeynCalc 10", arXiv:2312.14089 — FeynCalc conventions

**When citing, use this format in code:**
```julia
# Ref: refs/papers/Denner1993.pdf, Eq. (4.18)
# "B₁(p², m₀², m₁²) = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)"
```

## Build & test
```bash
# v2 tests (run individually — separate Julia processes)
for f in test/v2/test_*.jl; do julia --project=. "$f"; done

# v1 tests (616 tests, frozen — should always pass)
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Issue tracking: Beads + Dolt
All issues are tracked with Beads (not TodoWrite, not markdown files, not ad-hoc tracking).
```bash
# First time setup (or after clone)
bd init --force --prefix feynfeld && bd backup restore

# Check beads version — must be >= 0.61, must use Dolt backend
bd --version
bd doctor

# Workflow
bd ready                     # find available work
bd create --title="..." --description="..." --type=task --priority=2
bd update <id> --claim       # claim work
bd close <id>                # mark complete
bd dolt push                 # push to Dolt remote

# Session end: always push
bd dolt push
git push
```

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

1. **GROUND TRUTH = PHYSICS.** Not pinned numbers. Not LLM memory. Physics only.

2. **ANTI-HALLUCINATION: CITE EVERYTHING.** Every physics formula in source AND
   test code must cite: (a) local file path, (b) equation number, (c) verbatim
   copy of the equation. No uncited formulas. STOP and ask if source not found.
   ```julia
   # Ref: refs/papers/Denner1993.pdf, Eq. (4.18)
   # "B₁ = [A₀(m₁²) - A₀(m₀²) - (p²+m₀²-m₁²)B₀] / (2p²)"
   B1 = (a0_m1 - a0_m0 - (p2 + m02 - m12) * b0) / (2 * p2)
   ```

3. **SKEPTICISM.** Verify subagent work twice.

4. **ALL BUGS ARE DEEP.** Do not underestimate. No bandaids. Full solutions only.

5. **WORKFLOW.** 3 subagents before core code change (research source + 2 solutions).

6. **REVIEW.** Rigorous reviewer agent after every core change.

7. **TESTING.** Targeted only, or full suite in background.

8. **JULIA IDIOMATIC ALL THE WAY.** Read the cheatsheet below before writing code.

9. **NO PARALLEL JULIA AGENTS.** Read-only research agents CAN run in parallel.

10. **RESEARCH IDIOMS FIRST.** Before every new layer, research Julia patterns.

11. **LOC LIMIT.** No source file exceeds ~200 lines.

12. **REPEAT RULES.** Repeat occasionally to maintain focus.

## MUnit translation protocol
For each FeynCalc MUnit test file in `refs/FeynCalc/Tests/`:
1. Read the `.test` file.
2. Translate each `Test[]` to a Julia `@test`, preserving math exactly.
3. Document source file and test ID in a comment.
4. Cite the textbook equation that validates the test (Rule 2).
5. If MUnit test and textbook disagree, the textbook wins (Rule 1).

---

## Julia Idiom Cheatsheet

**READ THIS BEFORE WRITING ANY CODE.** Also in `JULIA_PATTERNS.md` and `Feynfeld_PRD.md` §6.

### DO

```julia
# Parametric types for dispatch (not tagged unions)
struct DiracGamma{S<:DiracSlot}; slot::S; end

# Multiple dispatch (not if/elseif isa cascades)
contract(a::MetricTensor, b::FourVector, idx) = ...
contract(a::FourVector, b::FourVector, idx) = ...

# Dict-based expression storage (O(1) like-term collection)
struct AlgSum; terms::Dict{FactorKey, Coeff}; end

# Concrete small Unions
const Coeff = Union{Rational{Int}, DimPoly}

# Named constructors (plain functions)
A0(m2::Real) = PaVe{1}(Int[], Float64[], [Float64(m2)])

# ScopedValues for implicit context
const CURRENT_SP = ScopedValue(SPContext())

# Use existing packages
using QuadGK: quadgk
using PolyLog: li2

# Dispatch, not type-in-function-name
evaluate(pv::PaVe{1}; mu2=1.0) = ...
evaluate(pv::PaVe{2}; mu2=1.0) = ...
```

### DO NOT

```julia
# ✗ isa checks instead of dispatch
if x isa ScalarProduct ... elseif x isa MetricTensor ... end

# ✗ Any/Expr coefficients
struct BadTerm; coeff::Any; end

# ✗ Type-in-function-name (Java)
evaluate_pave(x)       # ✗ → evaluate(x::PaVe)

# ✗ Wrapper structs for plain functions
struct Problem; alpha::Float64; end  # ✗ → f(; alpha=...)

# ✗ Hand-rolled numerics
const GAUSS_NODES = (...)  # ✗ → quadgk(f, 0, 1)

# ✗ Reimplementing standard functions
function my_dilog(z)...end  # ✗ → using PolyLog: li2

# ✗ Forcing types into wrong unions
AlgFactor = Union{..., PaVe}  # ✗ PaVe is scalar, not tensor

# ✗ Global mutable state
CURRENT = nothing  # ✗ → const CURRENT = ScopedValue(nothing)
```
