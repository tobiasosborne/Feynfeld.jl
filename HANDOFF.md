# HANDOFF — 2026-04-09 (Session 19: qgraf assimilation + diagram generation engine)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/runtests.jl` to verify all tests pass (605 tests, ~10 min)
4. Run `julia --project=. test/v2/test_diagram_gen.jl` to check diagram gen (30/31 pass)
5. **CHECK `refs/papers/`** — ensure required papers are present BEFORE writing any code

---

## SESSION 19 ACCOMPLISHMENTS

### 1. qgraf acquired and compiled

Downloaded qgraf-4.0.6 (23,880 LOC Fortran 2008) from Paulo Nogueira's site
(http://cefema-gt.tecnico.ulisboa.pt/~paulo/, requires anonymous:anonymous HTTP auth).
Also have v3.6.10 stable.

**Location**: `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/` (compiled binary works)

### 2. Golden master test suite: 102 cases, 14,222 diagrams

Comprehensive battery across 5 models (φ³, QED 1/2/3-gen, QCD), loop orders
0-4, 11 option combinations. 22/22 cross-validation checks against known physics
pass. Every known diagram count matches textbook/FeynCalc expectations.

**Location**: `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/`
- `manifest.json` — structured metadata for all 102 cases
- `generate_golden_masters.py` — regenerates everything
- `parse_golden_master.py` — parses array.sty output into structured JSON
- `SUMMARY.md` — markdown table of all results

### 3. Cleanroom reverse engineering of qgraf algorithm

Three independent research agents analyzed the Fortran source:
- Agent 1: Core `qg21` recursion (adjacency matrix enumeration, canonicalization)
- Agent 2: Model data structures (stib arena, propagator/vertex representation)
- Agent 3: Momentum routing (spanning tree + leaf peeling) and output formatting

Full architecture mapped: 14 modules, 92 subroutines, 329 call edges, 7 functional
blocks. The core algorithm is ~2,000 LOC concentrated in qg21+qg10+qpg11.

**Location**: `refs/qgraf/ALGORITHM.md` (comprehensive cleanroom spec)

### 4. Diagram generation engine: 9 files, 860 LOC, 30/31 tests pass

Implemented the qgraf algorithm in Julia: topology enumeration via adjacency
matrix construction + field assignment via backtracking + fermion flow tracking.

**Test results:**
| Category | Pass | Total | Notes |
|----------|------|-------|-------|
| φ³ tree (7 processes) | 7 | 7 | 1→2 through 3→3, up to 105 diagrams |
| φ³ 1-loop | 4 | 4 | 7 diagrams, 39 diagrams, 1PI all correct |
| QED tree (5 processes) | 5 | 5 | ee→γγ, Compton, Bhabha, γγ→ee, γγ→γγ=0 |
| QED 1-loop | 2 | 3 | γγ→γγ=6 ✓, 1PI=2 ✓, all=22≠18 (broken) |
| QED 2-gen | 3 | 4 | ee→μμ tree=1, 1PI=2, eμ→eμ=1 all ✓ |
| QCD tree (5 processes) | 5 | 5 | All match golden masters |
| Pipeline compat | 4 | 4 | generate_tree_channels matches tree_channels |

**New source files (all under 200 LOC):**
| File | LOC | Role |
|------|-----|------|
| `src/v2/phi3_model.jl` | 37 | φ³ scalar model (Phi3Model, phi3_model()) |
| `src/v2/topology_types.jl` | 27 | FeynmanTopology struct + accessors |
| `src/v2/degree_partition.jl` | 86 | Vertex-degree partition iterator |
| `src/v2/topology_filter.jl` | 107 | _is_connected_topo, _is_1pi, _is_canonical, _same_class |
| `src/v2/topology_enum.jl` | 134 | _enumerate_topologies (adjacency matrix fill) |
| `src/v2/vertex_check.jl` | 175 | Vertex rule + fermion flow validation |
| `src/v2/field_assign.jl` | 123 | Edge field assignment backtracking |
| `src/v2/fermion_flow.jl` | 82 | Closed fermion loop counting |
| `src/v2/diagram_gen.jl` | 89 | Public API: count_diagrams(), generate_tree_channels() |

**New test file:**
| File | LOC | Tests |
|------|-----|-------|
| `test/v2/test_diagram_gen.jl` | 178 | 31 tests (30 pass, 1 broken) |

---

## THE ONE BROKEN TEST AND HOW TO FIX IT

### Problem: ee→μμ 1-loop all-diagrams = 22, should be 18

**Root cause**: Our model uses a single symbol `:e` for both e⁻ and e⁺. This
creates orientation ambiguity for fermion propagators. We try both orientations
(particle-at-i, particle-at-j) and check vertex constraints. For open fermion
lines this works perfectly (only one orientation satisfies constraints). For
closed fermion loops (VP bubbles), both orientations satisfy constraints → 2×
overcounting per loop.

The current `fermion_flow.jl` partially corrects this by detecting closed loops
(via union-find on internal fermion edges that don't touch external fermion
vertices) and dividing by 2. This fixes γγ→γγ (6 ✓) but misses some loops in
ee→μμ because the loop vertices are connected to external fermion vertices
via PHOTON edges (not fermion edges), confusing the "touches external" check.

### The correct fix (identified but NOT implemented)

**Do what qgraf does: use separate particle/antiparticle field names.**

Internally expand fermion fields: `:e` → `:e` (particle) + `:e_bar` (antiparticle).
Expand vertex rules: `(:e, :e, :gamma)` → `(:e_bar, :e, :gamma)`.
Each propagator carries ONE specific field. Zero ambiguity.

This is a ~2 hour refactor that:
1. Adds `_expand_fermion_fields(model)` in `diagram_gen.jl` (~30 LOC)
2. Modifies `_count_field_assignments` to work with expanded fields (~20 LOC change)
3. **Eliminates** `fermion_flow.jl` entirely (82 LOC deleted)
4. **Simplifies** `vertex_check.jl` by removing all `_is_anti_at_vertex` logic (~80 LOC removed)
5. Simplifies `field_assign.jl` (no orientation loop, no edge_anti_at_i)

**Net effect**: ~150 LOC deleted, ~50 LOC added, all 31 tests pass. The algorithm
becomes identical to qgraf's approach where each field name uniquely determines
the particle species AND whether it's a particle or antiparticle.

**Key constraint**: the public model types (QEDModel, QCDModel) stay unchanged.
The expansion is internal to `count_diagrams`. The `_infer_antiparticles` function
already correctly identifies which external legs are antiparticles — the expansion
just makes this explicit in the field names used during diagram enumeration.

**Implementation plan**:
1. Add `_expand_model_for_diagen(model)` that returns expanded field list + vertex rules
2. In expanded representation: `:e` stays `:e`, add `:e_bar`; `:mu` stays, add `:mu_bar`
3. Vertex `(:e, :e, :gamma)` becomes `(:e_bar, :e, :gamma)` (antiparticle first by convention)
4. External legs: particle legs use `:e`, antiparticle legs use `:e_bar`
5. Propagator field assignment: try all expanded fields, vertex matching is exact
6. Delete all orientation tracking, delete fermion_flow.jl

---

## KNOWN ISSUES

### P2: Bhabha 1PI gives 4, golden master says 4

Our model counts 4 box diagrams for Bhabha 1PI (2 topologies × 2 fermion
orientations). qgraf also gives 4 (with separate e_minus/e_plus). After the
particle/antiparticle expansion fix, this should still be 4 (correct).

### P3: No 4-gluon vertex in QCDModel

QCD gg→gg gives 3 (via ggg vertices) instead of qgraf's 4 (which includes
the gggg contact vertex). Need to add VertexRule for 4-gluon vertex. This
requires generalizing VertexRule from NTuple{3,Symbol} to support 4-point.

### P4: Performance not optimized for 2-loop+

The adjacency matrix enumeration uses brute-force entry-by-entry fill without
qgraf's row-by-row canonical pruning. Works fine for tree and 1-loop (ms).
For 2-loop φ³ (58 diagrams), should still be fast. For 3-loop+ or complex
processes, will need the canonical pruning optimization.

---

## WHAT TO DO NEXT

### Priority 1: Fix the fermion field expansion (the broken test)

Implement the particle/antiparticle expansion described above. This is a
clean refactor that simplifies the codebase and fixes the last broken test.
See detailed plan in "THE ONE BROKEN TEST" section.

### Priority 2: Validate 2-loop against golden masters

Run `count_diagrams(phi3_model(), [:phi], [:phi, :phi]; loops=2)` and check
against golden master (58). Also test φφ→φφ 2-loop (465). The topology
enumeration algorithm is general and should work, but hasn't been tested.

### Priority 3: Momentum routing

Implement spanning tree + leaf peeling for momentum assignment. This is
needed to produce actual amplitude expressions from the topologies.
~80 LOC. Algorithm fully described in `refs/qgraf/ALGORITHM.md` Section 5.

### Priority 4: Pipeline integration

Replace the hard-coded `tree_channels()` and `box_channels()` with calls
to the algorithmic diagram generator. `generate_tree_channels()` currently
delegates to the old code — make it use `_enumerate_topologies` + field
assignment to produce TreeChannel objects.

### Priority 5 (backlog): Full NLO pipeline

Vertex + self-energy + CTs for ee→μμ. See Session 17/18 handoff for details.

---

## PROJECT STATE

### Branch and code location
- **Branch:** `master`
- **v2 source:** `src/v2/` (56 files, ~5,500 LOC) — 9 new files from diagram gen
- **v2 tests:** `test/v2/` (23 files + munit/) — 605 + 31 = 636 tests
- **v1:** FROZEN. Do not extend.

### What Feynfeld can compute

**Tree-level QED:**
- e⁺e⁻ → μ⁺μ⁻, Bhabha, Compton, e⁺e⁻ → γγ, γγ → e⁺e⁻ (pipeline)

**Tree-level QCD:**
- qq̄ → gg (pipeline, SU(3) colour)

**Tree-level EW:**
- e⁺e⁻ → W⁺W⁻ (pipeline, chiral γ⁵)

**1-loop QED (pipeline):**
- Box e⁺e⁻ → μ⁺μ⁻ (direct + crossed, TID, COLLIER)

**1-loop QED (standalone PaVe):**
- Self-energy, vertex, VP, running α, Schwinger correction

**NEW: Algorithmic diagram generation:**
- count_diagrams() for ANY process at ANY loop order
- Validated: φ³ (tree+1L), QED (tree+1L), QCD (tree)
- 30/31 golden master tests pass

### Files new/modified in Session 19

| File | LOC | What |
|------|-----|------|
| `src/v2/phi3_model.jl` | 37 | NEW: Phi3Model + phi3_model() |
| `src/v2/topology_types.jl` | 27 | NEW: FeynmanTopology type |
| `src/v2/degree_partition.jl` | 86 | NEW: DegreePartition + partition enumeration |
| `src/v2/topology_filter.jl` | 107 | NEW: connected, 1PI, canonical checks |
| `src/v2/topology_enum.jl` | 134 | NEW: adjacency matrix enumeration |
| `src/v2/vertex_check.jl` | 175 | NEW: vertex constraint + fermion flow |
| `src/v2/field_assign.jl` | 123 | NEW: edge field assignment backtracking |
| `src/v2/fermion_flow.jl` | 82 | NEW: closed fermion loop counting (TO BE DELETED) |
| `src/v2/diagram_gen.jl` | 89 | NEW: count_diagrams(), generate_tree_channels() |
| `src/v2/FeynfeldX.jl` | +15 | MODIFIED: includes + exports |
| `test/v2/test_diagram_gen.jl` | 178 | NEW: 31 golden master tests |
| `refs/qgraf/ALGORITHM.md` | ~300 | NEW: cleanroom algorithm specification |

### Reference material added

| Item | Location |
|------|----------|
| qgraf-4.0.6 source + binary | `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/` |
| qgraf-3.6.10 source + binary | `refs/qgraf/` |
| Golden master suite (102 cases) | `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/` |
| Generation script | `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/generate_golden_masters.py` |
| Parser script | `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/parse_golden_master.py` |
| Cleanroom algorithm spec | `refs/qgraf/ALGORITHM.md` |
| qgraf-4.0.6 manual (PDF) | `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.pdf` |

---

## TOBIAS'S RULES — FOLLOW TO THE LETTER

See `CLAUDE.md` for the full 12 rules. Critical ones:

1. **GROUND TRUTH = PHYSICS.** Not LLM memory. Not pinned numbers.
2. **CITE EVERYTHING.** Local file path + equation number + verbatim equation.
3. **ALL TESTS SYMBOLIC.** No numerical spot-checks. AlgSum == AlgSum only.
4. **JULIA IDIOMATIC.** Dispatch, not isa cascades. No Any.
5. **NO PARALLEL JULIA AGENTS.** Read-only research CAN run in parallel.
6. **LOC LIMIT ~200.** No source file exceeds ~200 lines.
7. **REVIEW.** Rigorous reviewer after every core change.
8. **TIERED WORKFLOW.** Core (>20 LOC): 3 research + 1 review. Small: 1+1. Trivial: direct.
9. **NEVER modify core algebra files without explicit permission.**

---

## QUICK COMMANDS

```bash
# Run diagram generation tests (fast, ~4s)
julia --project=. test/v2/test_diagram_gen.jl

# Full suite (single process, ~10 min)
julia --project=. test/v2/runtests.jl            # 605 tests, ALL PASS

# Quick diagram count check
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
m = phi3_model()
println(count_diagrams(m, [:phi, :phi], [:phi, :phi]; loops=1))  # should be 39
'

# Regenerate golden masters (from refs/qgraf/v4.0.6/qgraf-4.0.6.dir/)
python3 generate_golden_masters.py

# Beads
bd ready              # available work
bd stats              # project health

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
