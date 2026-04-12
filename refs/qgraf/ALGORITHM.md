# QGRAF 4.0.6 — Cleanroom Algorithm Specification

**Source**: qgraf-4.0.6.f08, 23,880 LOC Fortran 2008 (Paulo Nogueira, CeFEMA/IST)
**Analysis date**: 2026-04-09
**Purpose**: Specification for Julia reimplementation in Feynfeld.jl Layer 3
**Method**: Automated static analysis (fparser2) + manual source reading + 3 independent
research agents cross-validating. NO knowledge from external documentation used for the
algorithm description — this is purely from the source code.

---

## 0. Executive Summary

QGRAF generates all connected Feynman diagrams for a given process at a given
loop order. The algorithm has three phases:

```
Phase 1: PARTITION — iterate vertex-degree sequences (qpg11)
Phase 2: CONSTRUCT — build all non-isomorphic graphs per partition (qg21)
Phase 3: OUTPUT   — assign momenta via spanning tree, format output (qompac)
```

The core insight: qgraf builds topologies by **filling an adjacency matrix
row by row** with integrated canonicalization (isomorph rejection during
construction, not after). This avoids the exponential blowup of
generate-then-filter approaches.

**Architecture**: 14 modules (754 LOC data declarations), 92 subroutines
(~22,350 LOC). ~40% is I/O/parsing/error handling not needed in Julia.
The core algorithm is ~2,000 LOC concentrated in qg21+qg10+qpg11.

---

## 1. Data Structures

### 1.1 The Arena Allocator

All data lives in a single integer buffer `stib` (up to 8,388,607 entries)
and a character buffer `stcb`. "Pointers" are integer offsets into these
buffers. The sentinel `eoa` marks end-of-array; `nap` marks null.

```
stib(base_pointer + index) → value
stcb(position : position + length - 1) → string
```

### 1.2 Particle (Field) Representation

Each particle index `i` (1-based, 1..nphi) has parallel arrays:

| Array | Meaning |
|-------|---------|
| `namep(0)+i` | Position in stcb of name string |
| `namel(0)+i` | Length of name string |
| `link(0)+i` | Index of antiparticle (self if self-conjugate) |
| `antiq(0)+i` | Statistics: 0=boson, nonzero=fermion |
| `tpc(0)+i` | Type code: 0=default, 1=notadpole, 2=internal, 5=external |

For propagator `[A, B, +/-]`:
- If A=B (self-conjugate, e.g. photon): one slot, `link(i)=i`
- If A≠B (charged pair): two slots, `link(i)=i+1`, `link(i+1)=i`

### 1.3 Vertex Representation

Each vertex v has:
- `vval(0)+v` = degree (valence)
- `vparto(0)+v` = offset to sorted particle list in stib

The **`dpntro` lookup table** is the key structure for field insertion:
for each degree d, `dpntro(0)+d` points to a table indexed by first-particle.
All vertex configs of degree d are stored as sorted particle tuples,
enabling efficient iteration during field assignment.

### 1.4 Graph Representation

The graph is an **adjacency matrix** plus auxiliary arrays:

```
n            = total vertices (external + internal)
nli          = total lines (edges)
nleg         = number of external legs
rhop1        = nleg + 1 (first internal vertex index)

gam(i,j)    = number of edges between vertices i and j
               gam(i,i) = number of self-loops at vertex i
vdeg(i)     = degree of vertex i (1 for external, k for internal)
vmap(i,k)   = k-th neighbor of vertex i (with multiplicity)
lmap(i,k)   = back-pointer: if vmap(i,k)=j, then vmap(j,lmap(i,k))=i
amap(i,k)   = edge label for k-th half-edge at vertex i
pmap(i,k)   = field/propagator type at k-th half-edge at vertex i

head(l), tail(l) = endpoints of line l
intree(l)   = 1 if edge l is in spanning tree, 0 if chord
flow(l,0)   = edge type classification
flow(l,1..nleg) = external momentum coefficients
flow(l,nleg+1..nleg+nloop) = loop momentum coefficients
```

### 1.5 Five-Level Nesting

The generation has 5 nested levels, each with a control flag `qco`:

```
C (loop count)      → qg05: iterate nloop = loopx(1)..loopx(2)
  P (partition)     → qpg11: iterate vertex-degree partitions ρ
    T (topology)    → qg21: enumerate non-isomorphic adjacency matrices
      E (ext perms) → qg10: permute external legs, build vmap/lmap
        D (fields)  → qgen: assign fields to edges via dpntro tables
```

---

## 2. Phase 1: Vertex-Degree Partition (qpg11)

### 2.1 The Partition ρ

```
ρ(-1) = number of external vertices (= nleg, each degree 1)
ρ(k)  = number of internal vertices of degree k, k = mrho..nrho
```

Constraints:
```
V = Σ_k ρ(k) + ρ(-1)           (total vertices)
P = (Σ_k k·ρ(k) + ρ(-1)) / 2  (total half-edges / 2 = lines)
L = P - V + 1 = nloop          (Euler formula)
```

Equivalently: `Σ_k (k-2)·ρ(k) = nleg + 2·(nloop - 1)`

### 2.2 Iteration

For each loop count, qpg11 iterates all valid partitions:
1. Check `ρ(k) ≤ nivd(k)` (model has enough vertex types of degree k)
2. Check propagator feasibility (enough propagator types to fill lines)
3. For each valid partition, call qg21

---

## 3. Phase 2: Topology Construction (qg21) — THE CORE ALGORITHM

### 3.1 Overview

qg21 builds **labeled multigraph topologies** with prescribed vertex-degree
sequence. It does NOT enumerate all graphs then filter — it builds valid
topologies incrementally with integrated isomorph rejection.

Vertices are ordered: external first (1..ρ(-1)), then internal sorted by
degree (rhop1..n). The construction has three nested steps:

```
Step A: Distribute external connections (xn array)
  Step B: Distribute self-loops (diagonal of xg)
    Step C: Fill off-diagonal adjacency matrix row by row
```

### 3.2 Step A: External Leg Distribution

The ρ(-1) external legs are distributed among internal vertices:

```
xc(1) = external-to-external edges (must be even, step by 2)
xc(k) = total external legs going to vertices of degree degr(k)
xn(i) = number of external legs at vertex i

Iterate: for each xc(1) = 0, 2, 4, ..., ρ(-1):
  Distribute remaining ρ(-1)-xc(1) legs among internal vertex groups
  Within each group: xn(i) in non-increasing order (canonicalization)
```

External-to-external pairs are placed as: xg(2i-1, 2i) = 1 for i=1..xc(1)/2.

### 3.3 Step B: Self-Loop Distribution

```
For dsum = 0, 1, ..., nloop:
  dta(i) = max self-loops at vertex i (respects degree budget + options)
  Distribute 2·dsum self-loop half-edges among internal vertices
  xg(i,i) in non-increasing order within equivalence classes
  
  Backtrack (labels 32/33/65): find rightmost vertex that can accept more
```

### 3.4 Step C: Off-Diagonal Matrix Fill (THE HEART)

After fixing external connections and self-loops, fill the upper triangle
of the adjacency matrix row by row:

```
ds(lin, i) = remaining degree of vertex i at row lin
           = vdeg(i) - xn(i) - xg(i,i) - Σ_{k<lin} xg(k,i)
lps(lin)   = remaining loop budget at row lin
bond       = min(structural_bound, lps(lin))

For lin = limin to n:
  Fill xg(lin, col) for col = n, n-1, ..., lin+1 (greedy, largest first)
  xg(lin, col) = min(remaining_degree, bond, ds(lin, col))
  
  If remaining > 0 after all cols: BACKTRACK to previous row
  
  Check: msum = Σ (xg(lin,col)-1) for entries > 1
  If msum ≥ lps(lin): backtrack (not enough loops left)
```

### 3.5 Canonicalization (Isomorph Rejection)

**Equivalence classes** (xset/uset arrays): vertices with same
(degree, external-connection-count, self-loop-count) are interchangeable.

**Within-row ordering**: columns in the same equivalence class must have
non-increasing entries.

**Cross-row ordering**: if vertices i and i+1 are equivalent, row i must
be lexicographically ≥ row i+1.

**Permutation check** (after matrix completion): iterate all permutations
π of vertices within equivalence classes. Check gam(π(i),π(j)) vs gam(i,j):
- If gam(π(i),π(j)) > gam(i,j) for any pair: current labeling is NOT
  canonical → reject (backtrack)
- If identical for ALL pairs: π is an automorphism → count it (ngsym += 1)

This ensures each topology is generated exactly once in canonical form.

### 3.6 Connectedness Check

BFS from vertex n using signed coloring (aa array):
- Mark vertex n as +1
- Neighbors get opposite sign (-1)
- If all vertices reached: connected ✓
- If not all reached: disconnected → reject

The signed coloring also detects bipartiteness (for `bipart` option):
if any edge connects same-sign vertices → not bipartite.

### 3.7 Topological Filters

Applied AFTER matrix completion, BEFORE field insertion:

| Filter | Function | What it checks |
|--------|----------|----------------|
| `onepi`/`nobridge` | `qumpi(1,...)` | Remove each edge; BFS; if disconnects → bridge → reject |
| `nosbridge` | `qumpi(2,...)` | Bridge where one component has no external legs |
| `notadpole` | `qumpi(3,...)` | Bridge where one component has 0 or all external legs |
| `onshell` | `qumpi(4,...)` | Bridge where one component has exactly 1 external leg |
| `nosnail` | `qumvi(1,...)` | Self-loop on an internal vertex → reject |
| `onevi` | `qumvi(2,...)` | Remove each vertex; BFS; if disconnects → cut vertex → reject |
| `onshellx` | `qumvi(3,...)` | Cut vertex creating component with 1 external leg |
| `noselfloop` | In qg21 | Any xg(i,i) > 0 → reject |
| `nodiloop` | In qg21 | Any xg(i,j) > 1 for i≠j → reject |
| `noparallel` | In qg21 | Any xg(i,j) > 1 or xg(i,i) > 3 → reject |
| `cycli` | `qcyc` | Chord edges must form single connected cycle-block |
| `nosigma` | `qgsig` | No self-energy insertions (2-point subdiagrams) |

Filters with `< 0` flag REQUIRE the property (e.g., `selfloop` requires self-loops).

---

## 4. Field Insertion (qgen)

After qg10 returns a topology with external-leg permutation, qgen assigns
concrete field types to all edges.

### 4.1 Algorithm

```
External legs: pmap(i, 1) = leg(ps1(i)) for i = 1..nleg
               pmap(neighbor, back_slot) = link(leg(ps1(i)))

For each internal vertex vv in vlis order (connectivity-sorted):
  Look up vertex config via dpntro(degree) + first_assigned_particle
  For each candidate vertex configuration:
    Check: already-assigned legs match?
    Check: self-loop legs are valid (field, anti-field) pairs?
    Assign remaining legs, set conjugates at other endpoints
    Check: fermion flow consistency, exclusion filters
    
    If all vertices assigned: DIAGRAM COMPLETE → compute symmetry, output
    Else: advance to next vertex
  
  If no config works: BACKTRACK to previous vertex
```

### 4.2 External-Leg Permutations (qg10)

After each topology from qg21, qg10 iterates over all permutations of
external legs respecting symmetries. Skip permutations violating ordering
constraints for legs connected to the same internal vertex.

---

## 5. Momentum Routing

### 5.1 Spanning Tree (in qg21, after topology accepted)

Build spanning tree using greedy algorithm preferring high-multiplicity edges:

```
Initialize: all vertices unlabeled, no tree edges
For each multiplicity level (highest first):
  For each edge group:
    If endpoints in different components: merge, add to tree
    If both unvisited: start new component
    If same component: skip (would create cycle)
External legs always in tree. Stop when V-1 tree edges collected.
```

### 5.2 Momentum Assignment (Leaf Peeling)

```
Initialize:
  flow(leg_i, i) = 1 for i = 1..nleg  (each external leg carries pᵢ)
  flow(chord_j, nleg+k) = 1           (each chord carries loop momentum kₖ)

Repeat until all tree edges eliminated:
  Find leaf vertex (degree 1 in remaining tree)
  Its unique tree edge gets: flow = -(sum of flows of all other edges at leaf)
  Remove leaf from tree
  
Normalize: if majority of external coefficients are negative, flip signs
```

### 5.3 Edge Classification

```
flow(l, 0) = type:
  +1 (regular bridge):     only external momenta (all loop columns zero)
  +2 (special bridge):     zero momentum (all columns zero, tadpole-type)
  -1 (regular non-bridge): carries loop momenta (chord)
  -2 (special non-bridge): self-loop carrying only loop momenta
```

---

## 6. Symmetry Factors

**S = S_local × S_nonlocal**, output as 1/S.

### 6.1 Local Symmetry (S_local)

For each vertex pair (i,j) with multiplicity gam(i,j):
- **i ≠ j**: count groups of identical propagator assignments → k! per group
- **i = j** (self-loops): k! for groups of identical self-loop pairs,
  plus 2^k for self-conjugate self-loops (each can be flipped)

### 6.2 Nonlocal Symmetry (S_nonlocal)

Count automorphisms: permutations π within equivalence classes such that
gam(π(i),π(j)) = gam(i,j) for all pairs AND field assignment preserved.

---

## 7. Fermion Sign (qdis)

Trace all fermion lines through the diagram. Build a permutation from
canonical fermion ordering. Count transpositions → sign = (-1)^count.

Canonical order: incoming fermions 1,2,...,inco followed by
outgoing anti-fermions in reverse order.

---

## 8. Output Style Template

Four sections: `<prologue>`, `<diagram>`, `<epilogue>`, `<exit>`.

Key loop constructs in `<diagram>`:
- `<in_loop>` / `<out_loop>` — external legs
- `<propagator_loop>` — internal propagators
- `<vertex_loop>` + `<ray_loop>` — vertices and their half-edges
- `<momentum_loop>` — individual momentum terms

Key placeholders:
- `<field>`, `<momentum>`, `<field_index>`, `<propagator_index>`, `<vertex_index>`
- `<sign>`, `<symmetry_factor>`, `<diagram_index>`

**Field index convention** (signed):
- Positive odd: internal propagator, "from" end
- Positive even: internal propagator, "to" end
- Negative odd: incoming external leg
- Negative even: outgoing external leg

---

## 9. Mapping to Feynfeld.jl

| qgraf concept | Julia type | Notes |
|---------------|-----------|-------|
| stib arena | Native structs | No arena needed in Julia |
| Particle | `Field{Fermion/Boson}` | Already exists in Feynfeld |
| Propagator | `PropagatorRule` | Already exists |
| Vertex config | `VertexRule` | Already exists |
| dpntro table | `Dict{Int, Vector{VertexRule}}` | Degree → sorted configs |
| Adjacency matrix | `Matrix{Int}` or `SMatrix` | For small graphs |
| Degree partition | Iterator with `iterate()` | Julia iterator protocol |
| Topology | NEW: `FeynmanTopology` | Wraps adjacency matrix + metadata |
| Field assignment | `assign_fields(topo, model)` | Dispatch on (Topology, Model) |
| Momentum routing | `route_momenta(diagram)` | Spanning tree + leaf peeling |
| 1PI check | `is_1pi(topo)` | Bridge detection predicate |
| Symmetry factor | `symmetry_factor(diagram)` | Automorphism counting |
| Output | `build_amplitude(diagram)` | Already exists, extend |

### Key Julia idioms:

1. **No arena allocator** — Julia GC handles memory
2. **Adjacency matrix** — `Matrix{Int8}` for small graphs, StaticArrays for ≤20 vertices
3. **Partition iterator** — lazy `Channel` or custom iterator type
4. **Canonical ordering** — sort vertices by (degree, ext_count, self_loops)
5. **Filters** — composable predicates: `filter(is_1pi ∘ is_connected, topologies)`
6. **Field assignment** — multiple dispatch on `(Topology, QEDModel)` vs `(Topology, QCDModel)`
7. **Permutation iteration** — `Combinatorics.jl` or hand-rolled within equivalence classes

### Estimated Julia LOC:

| Component | LOC | Notes |
|-----------|-----|-------|
| Partition iterator | ~60 | Replaces qpg11 (248 LOC) |
| Adjacency matrix builder | ~200 | Replaces qg21 core (1242 LOC) |
| Canonicalization | ~80 | Equivalence classes + permutation check |
| Topological filters | ~100 | is_1pi, is_connected, has_tadpole, etc. |
| Field assignment | ~100 | Vertex config matching |
| Momentum routing | ~80 | Spanning tree + leaf peeling |
| Symmetry factor | ~40 | Automorphism counting |
| **Total** | **~660** | vs ~8,000 LOC in Fortran (excluding I/O) |

---

## 10. Golden Master Test Suite

102 cases, 14,222 diagrams, 22/22 cross-validations pass.

Located in: `refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/`

| Model | Cases | Diagrams | Max loops | Largest case |
|-------|-------|----------|-----------|-------------|
| φ³ | 40 | 13,753 | 4 | φφ→φφ 3L (5,625) |
| QED 1-gen | 29 | 271 | 2 | e→eγ 2L (50) |
| QED 2-gen | 11 | 44 | 1 | e+e-→μ+μ- 1L (18) |
| QED 3-gen | 6 | 39 | 2 | e+e-→μ+μ- 1L (23) |
| QCD | 16 | 115 | 1 | gg→ggg 0L (25) |

Key validation points:
- e+e-→μ+μ- tree: 1 diagram (s-channel) ✓
- e+e-→μ+μ- 1L 1PI: 2 diagrams (direct + crossed box) ✓
- e+e-→μ+μ- 1L all: 18 diagrams (2gen) / 23 diagrams (3gen) ✓
- γγ→γγ tree: 0 diagrams (Furry's theorem) ✓
- γγ→γγ 1-loop: 6 diagrams (box permutations) ✓
- Bhabha 1PI 1L: 4 diagrams (2 channels × 2 boxes) ✓
- γ self-energy 1PI (3gen): 3 diagrams (e,μ,τ bubbles) ✓
