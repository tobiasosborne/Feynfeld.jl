# HANDOFF — 2026-04-10 (Session 20: fermion field expansion fix + topology audit)

## DO NOT DELETE THIS FILE. Read it completely before working.

---

## START HERE

1. Read `CLAUDE.md` — rules, **pipeline principle**, anti-hallucination, Julia idioms
2. Run `bd ready` to see available work
3. Run `julia --project=. test/v2/test_diagram_gen.jl` to verify diagram gen (31/31 pass)
4. Run `julia --project=. test/v2/runtests.jl` to verify full suite (605 tests, ~5 min)
5. **READ THE CRITICAL FAILURE SECTION BELOW BEFORE DOING ANYTHING**

---

## SESSION 20 ACCOMPLISHMENTS

### 1. Fermion field expansion fix (the broken test)

Fixed `count_diagrams(qed, [:e,:e], [:mu,:mu]; loops=1)` returning 22 instead
of 18. Replaced the old orientation-tracking approach (`edge_anti_at_i` +
`_check_vertex_with_flow` + `fermion_flow.jl`) with qgraf's expanded-field
approach: separate `:e`/`:e_bar` symbols for particle/antiparticle.

**Key insight the Session 19 agent missed and the stuck Session 20a agent got
wrong**: when an edge carries field `f`, vertex i (lower index) sees `f` but
vertex j (higher index) sees `conjugate(f)`. The stuck agent used `f` at both
ends — that's why it broke 10 tests.

**Two overcounting corrections needed with expanded fields:**
1. **Parallel edge canonical ordering**: multi-edges between same vertex pair
   must have fields in non-decreasing order (prevents permutation overcounting)
2. **Closed fermion loop ÷2**: self-loops (edge from v to v) with fermion
   fields are counted in both orientations; skip multi-edge components
   (already handled by #1)

**Files changed:**
| File | Before | After | Change |
|------|--------|-------|--------|
| `diagram_gen.jl` | 89 LOC | 149 LOC | +ExpandedModel, expansion functions, -_infer_antiparticles |
| `field_assign.jl` | 123 LOC | 152 LOC | Rewritten: expanded backtracking + loop correction |
| `vertex_check.jl` | 175 LOC | 74 LOC | Rewritten: conjugate-aware vertex checks, deleted flow functions |
| `fermion_flow.jl` | 82 LOC | DELETED | Replaced by _count_closed_loops_expanded |
| `FeynfeldX.jl` | — | -1 line | Removed fermion_flow.jl include |
| `test_diagram_gen.jl` | — | -2 lines | @test_broken → @test |

**Test results: 31/31 diagram gen, 605/605 full suite.**

### 2. Comprehensive golden master validation

Ran 46 golden master cases against qgraf counts (up to 2-loop):
- **44 PASS** — all tree-level, all 1-loop, all QED, all QCD
- **2 KNOWN FAIL** — 2-loop φ³ (474 vs 465): duplicate topologies from
  incomplete canonicalization (see CRITICAL FAILURE below)

### 3. Gap analysis and issue registration

Registered 7 beads issues for golden master coverage gaps:
- `feynfeld-ejl` P2: 1-gen QED model (unlocks 13 cases)
- `feynfeld-wjw` P2: 4-point VertexRule for gggg (unlocks 8 cases)
- `feynfeld-rff` P3: 3-gen QED model (unlocks 6 cases)
- `feynfeld-sia` P3: diagram filter predicates (unlocks 13 cases)
- `feynfeld-jzg` P4: QCD ghost fields (unlocks 1 case)
- `feynfeld-ciz` P2: **topology_enum.jl rewrite** (correctness + performance)
- `feynfeld-0b2` P2: 2-loop overcounting (symptom of feynfeld-ciz)

---

## !! CRITICAL FAILURE: topology_enum.jl IS A BROKEN PROTOTYPE !!

**This is the most important section of this handoff.**

`topology_enum.jl` does NOT implement qgraf's algorithm. It implements a naive
brute-force approximation that happens to give correct results at tree and
1-loop but **produces wrong answers at 2-loop+**.

### What qgraf does (ALGORITHM.md §3):

```
Step A: Distribute external connections → xn(i), non-increasing within equiv classes
Step B: Distribute self-loops → xg(i,i), non-increasing within equiv classes
Step C: Fill off-diagonal ROW BY ROW, largest column first
        → within-row: non-increasing within equivalence class
        → cross-row: row i ≥ row i+1 lexicographically for equivalent vertices
        → loop budget checked at each row
        → full permutation-class isomorph check after completion
```

Each step prunes entire search tree branches BEFORE entering the next level.

### What our code does:

```
For each entry (i,j) in upper triangle (flat iteration):
    Try k = max_val, max_val-1, ..., 0
    Recurse to next entry
When ALL entries filled:
    Check connected? canonical? 1PI?
    Keep if all pass
```

**Three fatal defects:**

1. **No row-level pruning** — generates exponentially many invalid matrices,
   rejects at the end. A 3-loop φ→φφ that qgraf handles in milliseconds
   takes minutes (stack depth >80 in `_fill_entry!`).

2. **Incomplete canonicalization** — `_is_canonical_topo` only checks pairwise
   vertex swaps. qgraf checks ALL permutations within equivalence classes
   (§3.5). At 2-loop+ with 3+ equivalent vertices, pairwise swaps miss
   isomorphisms → **duplicate topologies → wrong diagram counts** (474 vs 465).

3. **No structural constraints during fill** — no loop budget tracking, no
   propagator feasibility, no row ordering. All checked post-hoc.

### What must be done (feynfeld-ciz):

Port qgraf's actual qg21 algorithm from ALGORITHM.md §3:
1. Step A: external leg distribution with equivalence-class ordering (~40 LOC)
2. Step B: self-loop distribution with ordering (~30 LOC)
3. Step C: row-by-row off-diagonal fill with within-row and cross-row
   ordering + loop budget (~80 LOC)
4. Full permutation-class isomorph check (~50 LOC)

Total: ~200 LOC Julia replacing the current 134 LOC. This is NOT an
optimization — it is a **correctness fix** for 2-loop+ topology enumeration.

**Do NOT attempt to patch the current approach.** The entry-by-entry fill
with post-hoc validation is fundamentally the wrong algorithm. Port qg21
properly or the 2-loop counts will never be right.

---

## WHAT TO DO NEXT

### Priority 1: Port qg21 properly (feynfeld-ciz + feynfeld-0b2)

This is the #1 task. Everything else is blocked by correct topology
enumeration. The algorithm is fully specified in ALGORITHM.md §3. The
golden master suite (102 cases, 14k diagrams) provides complete validation.
The current 31 tests in test_diagram_gen.jl MUST continue to pass (they're
all tree+1-loop which the current code handles correctly).

**Approach:**
1. Read ALGORITHM.md §3 completely (Steps A, B, C, canonicalization)
2. Write new `topology_enum_v2.jl` implementing qg21 properly
3. Red-green TDD: start with the 2-loop cases that currently fail
4. Validate against ALL golden masters
5. Replace `topology_enum.jl` once all tests pass

### Priority 2: Close golden master coverage gaps

After topology enumeration is correct:
- Add 1-gen QED model (feynfeld-ejl, ~10 LOC, unlocks 13 cases)
- Generalize VertexRule to N-point (feynfeld-wjw, ~30 LOC, unlocks 8 cases)
- Add 3-gen QED model (feynfeld-rff, ~5 LOC, unlocks 6 cases)
- Add diagram filter predicates (feynfeld-sia, ~100 LOC, unlocks 13 cases)

### Priority 3: Momentum routing

Implement spanning tree + leaf peeling for momentum assignment.
~80 LOC. Algorithm in ALGORITHM.md §5. Needed for actual amplitudes.

### Priority 4: Pipeline integration

Replace hard-coded tree_channels() with algorithmic generation.

---

## WHAT WORKS (do not break)

| Component | Status | Tests |
|-----------|--------|-------|
| Expanded field assignment | CORRECT at all loop orders | 31/31 |
| φ³ tree + 1-loop | CORRECT | 11/11 golden masters |
| QED tree (all processes) | CORRECT | 14/14 golden masters |
| QED 1-loop 1PI | CORRECT | 7/7 golden masters |
| QED 1-loop all (2-gen) | CORRECT (ee→μμ=18 ✓) | 2/2 golden masters |
| QCD tree (3-point only) | CORRECT | 5/5 golden masters |
| Full v2 algebra suite | CORRECT | 605/605 |

## WHAT IS BROKEN

| Component | Status | Issue |
|-----------|--------|-------|
| 2-loop+ topology enumeration | WRONG COUNTS | feynfeld-ciz |
| 2-loop+ canonicalization | INCOMPLETE (pairwise only) | feynfeld-0b2 |
| 3-loop+ performance | UNUSABLE (minutes+) | feynfeld-ciz |
| QCD 4-gluon vertex | MISSING | feynfeld-wjw |
| 1-gen QED model | MISSING | feynfeld-ejl |
| Diagram filter options | MISSING | feynfeld-sia |

---

## PROJECT STATE

### Branch and code location
- **Branch:** `master`
- **v2 source:** `src/v2/` (55 files, ~5,400 LOC) — fermion_flow.jl deleted
- **v2 tests:** `test/v2/` (23 files + munit/) — 605 + 31 = 636 tests
- **v1:** FROZEN. Do not extend.

### What Feynfeld can compute (unchanged from Session 19)

**Tree-level QED/QCD/EW:** e⁺e⁻→μ⁺μ⁻, Bhabha, Compton, e⁺e⁻→γγ, qq̄→gg, e⁺e⁻→W⁺W⁻
**1-loop QED (pipeline):** Box e⁺e⁻→μ⁺μ⁻ (direct + crossed, TID, COLLIER)
**1-loop QED (standalone):** Self-energy, vertex, VP, running α, Schwinger
**Algorithmic diagram generation:** count_diagrams() for ANY process, CORRECT at tree+1-loop

### Files changed in Session 20

| File | LOC | What |
|------|-----|------|
| `src/v2/diagram_gen.jl` | 149 | MODIFIED: +ExpandedModel, expansion functions |
| `src/v2/field_assign.jl` | 152 | REWRITTEN: expanded backtracking + loop correction |
| `src/v2/vertex_check.jl` | 74 | REWRITTEN: conjugate-aware, old flow functions deleted |
| `src/v2/fermion_flow.jl` | — | DELETED (82 LOC) |
| `src/v2/FeynfeldX.jl` | -1 | Removed fermion_flow.jl include |
| `test/v2/test_diagram_gen.jl` | -2 | @test_broken → @test |

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
8. **TIERED WORKFLOW.** Core (>20 LOC): 3 research + 1 review.

---

## QUICK COMMANDS

```bash
# Diagram generation tests (fast, ~3s)
julia --project=. test/v2/test_diagram_gen.jl

# Full suite (single process, ~5 min)
julia --project=. test/v2/runtests.jl

# Quick diagram count check
julia --project=. -e '
include("src/v2/FeynfeldX.jl"); using .FeynfeldX
println(count_diagrams(qed_model(), [:e, :e], [:mu, :mu]; loops=1))  # should be 18
'

# Beads
bd ready              # available work
bd stats              # project health
bd show feynfeld-ciz  # THE critical issue

# Session end protocol
git add <files> && git commit -m "..." && git push
bd dolt push
```
