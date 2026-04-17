# STOCKTAKE: Feynfeld.jl qgraf Port (Session 27, Phase 18b kickoff)

**Date:** 2026-04-17  
**Scope:** `src/v2/qgraf/` — 15 files, 3583 LOC  
**Context:** Strategy C qgraf-4.0.6 port; Phase 17 (audition/canonical) complete; Phase 18a (tree amplitudes) shipped; Phase 18b-1 skeleton landed (Session 27).

---

## Per-File Summary Table

| File | LOC | Purpose | Exports | qgraf refs | Phase 18b Notes |
|------|-----|---------|---------|------------|-----------------|
| **QgrafPort.jl** | 61 | Module entry; includes & exports all | 27 public: Partition, EquivClass, FilterSet, TopoState, step_b_enumerate!, step_c_enumerate!, qg21_enumerate!, etc. | — | Export list extended: added `combine_m_squared_burnside`, `is_emission_canonical`, dedup counters |
| **types.jl** | 194 | Core data structures: Partition, EquivClass, FilterSet, TopoState | Partition, EquivClass, FilterSet, TopoState, rho_k, n_internal, n_vertices, n_edges, no_filters, MAX_V | qgraf-4.0.6.f08:1743–1991 (qpg11 Partition) | TopoState defers uset, str, xl, xt, xp, etc. to "later phases"; Phase 18a-3 fills them in |
| **topology.jl** | 744 | **Biggest file:** Step B (xc/xn distribution) + Step C (xg topology generation) + helpers | step_b_enumerate!, _degree_class_bounds, _is_connected_internal, qg21_enumerate!, qg10_enumerate!, step_c_enumerate! | qgraf-4.0.6.f08:12426–13150 (qg21 full) | Core recursive-descent engine for topology generation; minimal phase 18b involvement (filters stay clean) |
| **canonical.jl** | 413 | Full per-equiv-class permutation canonicalisation (Knuth Algorithm L per class) | compute_equiv_classes!, next_class_perm!, is_canonical_full!, is_canonical_qgraf!, is_canonical_feynman, enumerate_topology_automorphisms, build_dpntro | qgraf-4.0.6.f08:12819–12840, 13156–13291 | Fixed 474→465 topology bug (pairwise-swap → full-permutation). Automorphism enumeration used by Phase 17 audition & Phase 18b Burnside weight calc |
| **filters.jl** | 234 | Topological filter predicates (selfloop, diloop, parallel, bridge, tadpole, etc.) | has_no_selfloop, has_no_diloop, has_no_parallel, is_one_pi, has_no_sbridge, has_no_tadpole, has_no_onshell, has_no_snail, is_one_vi | qgraf-4.0.6.f08:3690–3776 (qumpi), 3777 (qumvi), 13669 (qgsig), 18830 (qcyc) | Phase 14 baseline filters. Ported 9/14 qgraf options. No phase 18b blocker. |
| **audition.jl** | 339 | Phase 17 dedup audit: three strategies (Burnside, canonical-pmap, pre-filter) | count_dedup_burnside, count_dedup_canonical, count_dedup_prefilter, is_emission_canonical, emission_stabilizer, pmap_stabilizer, _foreach_emission, count_diagrams_qg21 | qgraf-4.0.6.f08:14036–14100 (qgen canonical-pmap dedup) | **KNOWN BUG:** canonical filter rejects one Bhabha orbit (vjw9 blocker); Burnside vs canonical disagree on ee→ee (expected, strategy-specific) |
| **qgen.jl** | 612 | Field assignment via dpntro lookup + recursive backtracker + slot-perm enumeration | build_dpntro, compute_qg10_labels, qgen_count_assignments, qgen_enumerate_assignments, qdis_fermion_sign, compute_local_sym_factor | qgraf-4.0.6.f08:13797–14464 (full qgen + filter suite) | Phase 12 complete. **TODO BUG 2:** pmap[vv,rdeg+1..vdeg] not saved on backtrack (φ³ 2-loop self-loop topologies over-count); mitigated by filters. |
| **spanning_tree.jl** | 69 | BFS spanning tree + chord count (Phase 16 minimal) | build_spanning_tree, count_chords | qgraf-4.0.6.f08:13315–13402 | Tree/chord classification. Chord count = nloop (Euler). No phase 18b work. |
| **momentum.jl** | 230 | Leaf-peel momentum routing (Phase 18a-1) | InternalEdge, EdgeMomenta, route_momenta | qgraf-4.0.6.f08:13400–13559 | Assigns momenta to internal edges under qgraf "all incoming" convention. No phase 18b scope. |
| **halfedge.jl** | 89 | Half-edge global labelling compute_amap (Phase 18a-2) | compute_amap | qgraf-4.0.6.f08:12133–12158, 12342–12344, 13320–13340 | Edge ids; pairing invariant. Used by all Phase 18a-3+ assemblers. No phase 18b work. |
| **propagator_assemble.jl** | 104 | Per-edge propagator factors (Phase 18a-3) | Propagator, build_propagators | qgraf-4.0.6.f08:— | Tree-level massless propagators (Boson: 1, Fermion: p̸+m, Scalar: 1). **Deferred:** composite/zero-momentum fermion propagators to 18b. |
| **vertex_assemble.jl** | 182 | Per-vertex Lorentz factors (Phase 18a-4) | build_vertices | qgraf-4.0.6.f08:— | 3-vertex QED/QCD/EW Lorentz index naming `:mu_l_<edge_id>` (mirrors src/v2/amplitude.jl convention). **Error on 4-vertex:** deferred to 18b. Index unification needed when explicit metrics appear. |
| **fermion_line.jl** | 93 | Fermion-line traversal (Phase 18a-6) | FermionLine, walk_fermion_lines | qgraf-4.0.6.f08:— | Tree-only: single vertex per line. **Error on internal fermion propagators:** deferred to 18b (Compton tree, fermion loops). |
| **emission_amplitude.jl** | 137 | Master assembler emission_to_amplitude (Phase 18a-7) | AmplitudeBundle, emission_to_amplitude | qgraf-4.0.6.f08:— | Composes phases 18a-1..6. Tree scope; bails on internal-fermion errors. Coupling placeholder alg(1) → phase 18b-7. |
| **burnside_combine.jl** | 82 | **New, Session 27:** Multi-orbit Burnside |M|² summation (Phase 18b-1) | combine_m_squared_burnside | — | **Option A scope:** trace-only AlgSum; 1/denom lives in bundle.denoms, caller applies 1/denom. Diagonal (spin_sum_amplitude_squared) + off-diagonal (spin_sum_interference). Fermion signs via qdis_fermion_sign. Supports 0- and 2-line bundles only; multi-vertex lines (18b-3). |

---

## Pipeline Position

The qgraf port implements **Layer 3 (Diagrams)** of the 12-layer Feynfeld pipeline:

1. **Input:** QED/QCD/EW model (Feynman rules, field species, vertex arities).  
2. **Step A:** Partition enumeration (Phase 10, external legacy).  
3. **Step B–C:** qg21 topology enumeration → canonical topologies (xg adjacency).  
4. **Step B–C scope:** 3583 LOC spread across types, topology recursion, canonicalisation, filters, qgen field assignment.  
5. **Step D:** Per-emission amplitude assembly (phases 18a-1..7):
   - route_momenta (leaf-peel algorithm).
   - compute_amap (half-edge labelling).
   - build_propagators (p/m² and numerators).
   - build_vertices (Lorentz factors, 3-vertex only).
   - walk_fermion_lines (tree line pairs).
   - build_externals (spinor/polarisation slots).
   - emission_to_amplitude (composition).
6. **Step E:** Multi-orbit combination via Burnside (Phase 18b-1, Option A trace-only).  
7. **Output:** AmplitudeBundle per canonical emission → Layer 4 (cross-section + spin-sum).

**Key property:** Zero allocations in hot path; all scratch arrays pre-allocated to MAX_V=24 at TopoState construction.

---

## Known Gaps & Phase 18b Deferrals

### Blocker: feynfeld-vjw9 (18b-1a, canonical orbit-rep dedup)

**Session 27 discovery:** `solve_tree_pipeline(Bhabha ee→ee)` reports `n_emissions=1` but qgraf correctly counts 2 (s + t orbits). The `is_emission_canonical` filter in `audition.jl:69` rejects one of the two orbits.

**Root cause:** Strategy C canonical-pmap check is valid for topology automorphisms but may be **invalid for qgen orbit-reps**. The orbit-rep may differ in pmap signature under qgen's flavor-assignment logic, causing a canonical-rep of one orbit to fail the filter. This is the bug documented in HANDOFF Session 22 Phase 17a VERDICT.

**Impact:** Phase 18a regression green (ee→μμ canonical-filter-compatible); Bhabha validation blocked.

**File:** `audition.jl:69–86` (`is_emission_canonical` function).

### Phase 18b-1 (trace-only AlgSum, Option A)

**Implemented (Session 27):** `burnside_combine.jl` with dual-trace formula:  
Σ_i,j w_i w_j s_i s_j T_ij (diagonal spin_sum_amplitude_squared, off-diagonal spin_sum_interference).

**Deferred (Option B):** Inverse-denominator factor `InverseSP` for symbolic 1/pair(q,q); tracked in bead feynfeld-rj1l.

### Phase 18b-2 (composite-momentum fermion propagators)

**File:** `propagator_assemble.jl:97`  
**Issue:** Fermion propagators with composite/zero momentum error.  
**Status:** Placeholder mass (0//1 or 1//1) only; symbolic mass support deferred.

### Phase 18b-3 (multi-vertex fermion lines)

**Files:** `fermion_line.jl:65`, `emission_amplitude.jl:51`, `burnside_combine.jl:73`  
**Issue:** Compton tree (2+ internal vertices per line) and fermion loops not supported.  
**Current:** Tree-only, single-vertex lines per QED 3-vertex.

### Phase 18b-4 (4-vertex, all-boson)

**File:** `vertex_assemble.jl:48`  
**Issue:** 4-vertex gggg (QCD, EW) errors with deferral message.  
**Status:** Placeholder only.

### Phase 18b-5 (boson polarisation)

**File:** `vertex_assemble.jl:114`, `emission_amplitude.jl:5`  
**Issue:** Explicit boson polarisation sums deferred; currently implicit via contraction.  
**Index naming:** `:mu_l_<edge_id>` (this file) vs `:mu_<channel>` (handbuilt amplitude.jl). Both collapse under spin-sum but **unification needed** when explicit metrics or multi-orbit interference require free boson indices.

### Phase 18b-7 (coupling assignment)

**File:** `emission_amplitude.jl:27`  
**Issue:** Coupling placeholder `alg(1)` only; full e², α_s, etc. deferred.

### Phase 18b-8 (validation)

**Test:** `test/v2/qgraf/test_phase18b1_multi_orbit.jl` (not yet written).  
**Gating:** Awaits vjw9 fix.

---

## Technical Debt & Known Bugs

| File | Issue | Severity | Workaround |
|------|-------|----------|-----------|
| `canonical.jl:3–15` | Old pairwise-swap canonicalisation (src/v2/topology_enum.jl) was bug locus; generated 474 duplicates at φ³ 2L instead of 465. | **Fixed in Phase 17.** | Full-permutation Knuth Algorithm L per equiv-class. |
| `qgen.jl:588–592` | **TODO BUG 2:** pmap[vv,rdeg+1..vdeg] not saved on backtrack (φ³ 2-loop self-loop topologies over-count). | **Medium.** Mitigated by filters; affects only self-loop-heavy topologies. | Existing filter suite (nosl, nodl, nopa) rejects most problematic cases. |
| `audition.jl:69` | **vjw9 blocker:** `is_emission_canonical` rejects valid Bhabha orbit (canonical orbit-rep invalid for qgen). | **High.** Blocks Bhabha 18b validation. | Awaits investigation: may need orbit-rep-aware canonical check or qgen re-canonicalisation. |
| `vertex_assemble.jl:67` | Boson Lorentz index naming mismatch (`:mu_l_<edge_id>` here vs `:mu_<channel>` in amplitude.jl). | **Medium.** Currently masked by spin-sum contraction. | Unification mandatory when explicit metrics or free-index interference appear in 18b. |
| `propagator_assemble.jl:50` | PaVe machinery & explicit `-g^{μν}` metric deferred to 18b. | **Expected phase boundary.** | Placeholder rational mass (1//1). |

---

## Code Quality Notes

- **Zero-allocation hot path:** All 15 files honor the pre-allocation discipline (MAX_V=24 scratch arrays in TopoState).
- **Refactoring coverage:** qgraf goto-style control flow faithfully mirrored in Julia @label/@goto (cross-auditable line-by-line).
- **Export discipline:** All 27 public functions listed in QgrafPort.jl module header.
- **Test coverage:** Phase 17 audition (audition.jl) passes; Phase 18a green (ee→μμ tree); Phase 18b-1 skeleton green except vjw9 (Bhabha).

---

## Session 27 Handoff Summary (150 words)

**Phase 18b kickoff landed with Option A (trace-only) Burnside combine and the vjw9 blocker (canonical orbit-rep dedup).** The skeleton (`burnside_combine.jl`, 82 LOC) composes `spin_sum_amplitude_squared` (diagonal) and `spin_sum_interference` (off-diagonal) under dual-Burnside weights and fermion signs. Phase 18a (tree amplitudes) stays green for ee→μμ; Bhabha fails because `is_emission_canonical` in `audition.jl:69` rejects one of two orbits. Root cause: the canonical-pmap check assumes topology automorphisms preserve qgen's orbit structure, but it may not—exactly the Strategy-C bug from HANDOFF Session 22. **Next agent:** start with feynfeld-vjw9, debug the canonical filter, then write `test_phase18b1_multi_orbit.jl` to gate Bhabha acceptance. Full Phase 18b roadmap spans 8 subtasks (18b-2 composite momenta, 18b-3 multi-vertex fermion lines, 18b-4 4-vertex, 18b-5 polarisation, 18b-7 coupling, 18b-8 validation) wired as epic `feynfeld-xa7s` with dependency graph.

