# Stocktake: v2 Model/Rules/Diagrams (Layers 1-3)

**Date:** 2026-04-17  
**Scope:** Model (Layer 1), Rules (Layer 2), Diagrams (Layer 3)  
**Total LOC (surveyed):** 1,307 lines  
**Status:** Layers 1-2 stable; Layer 3 in active development (Phase 18a-7)

---

## File-by-File Summary

| File | LOC | Public Exports | Purpose |
|------|-----|-----------------|---------|
| **Layer 1: Model** |
| `model.jl` | 116 | `GaugeGroup`, `Field{S}`, `AbstractModel`, `QEDModel`, `qed_model()`, `qed1_model()`, `qed3_model()`, `FieldSpecies`, `Fermion`, `Boson`, `Scalar` | Abstract gauge-group & field interface; concrete QED (e+e⁻→γ) model with 1-3 generations |
| `ew_model.jl` | 49 | `EWModel`, `ew_model()`, `feynman_rules(::EWModel)` | Standard Model electroweak for tree-level e⁺e⁻→W⁺W⁻ (5 vertices: eeγ, eeZ, eνW, WWγ, WWZ) |
| `ew_parameters.jl` | 39 | `EW_M_W`, `EW_M_Z`, `EW_ALPHA`, `EW_SIN_W`, `EW_GV_E`, `EW_GL_E` (+ Rational versions) | PDG 2024 on-shell scheme: masses, couplings, derived ratios for ew_model |
| `qcd_model.jl` | 81 | `QCDModel`, `qcd_model()`, `feynman_rules(::QCDModel)`, `triple_gauge_vertex()` | Quark-gluon + 3/4-gluon vertices; ghost field (modeled as fermion); momentum-dependent triple-gluon Lorentz structure |
| `phi3_model.jl` | 37 | `Phi3Model`, `phi3_model()`, `feynman_rules(::Phi3Model)` | Scalar φ³ for pure topology testing (no fermion flow, no color) |
| **Layer 2: Rules** |
| `rules.jl` | 121 | `VertexRule`, `FeynmanRules`, `vertex_factor()`, `propagator_num()`, `vertex_structure()`, `gauge_coupling_phase()`, `arity()` | Callable struct dispatch on field species & coupling keys; Lorentz structures (γ^μ, γ^μ(gV−gAγ5), γ^μ(1−γ5)/2); propagator numerators; generic fallback rule generator |
| **Layer 3: Diagrams** |
| `diagrams.jl` | 22 | `ExternalLeg`, constructor + show | Particle momentum/spin state at 4-leg interface (backward-compat constructor: mass→0) |
| `diagram_gen.jl` | 177 | `count_diagrams()`, `generate_tree_channels()`, `_count_diagrams_legacy()` | **Phase 17c:** delegates to `QgrafPort.count_diagrams_qg21`; legacy in-Julia enum preserved for regression; expands fermion fields (particle≠antiparticle) |
| `topology_enum.jl` | 117 | `_enumerate_topologies()`, `_fill_entry!()` | Adjacency-matrix backtracking (upper-triangle fill) with dedup via `QgrafPort.is_canonical_feynman` (qgraf lex-next permutation, qg21 labels 77/93/102/202/204) |
| `topology_filter.jl` | 107 | `_is_connected_topo()`, `_is_1pi()`, `_is_canonical()`, `_same_class()` | Connectedness check (BFS); 1PI validation (internal bridge detection); vertex-equivalence class dedup (pairwise swap lex-compare) |
| `topology_types.jl` | 27 | `FeynmanTopology`, `n_vertices()`, `n_internal()`, `n_propagators()`, `n_loops()` | Minimal type: n_ext, adjacency matrix (Int8), vertex degrees; Euler formula: L = P−V+1 |
| `channels.jl` | 103 | `TreeChannel`, `tree_channels()`, `vertex_legs()` | Tree 2→2 s/t/u hard-coded (filtering-only, no qgraf needed); validates vertex existence + fermion-pair spinor ordering (bar left, plain right) |
| `degree_partition.jl` | 86 | `DegreePartition`, `_degree_partitions()`, `_model_vertex_degrees()` | Recursive partition enumeration: Σ(k−2)ρ(k) = n_ext + 2(L−1); validates Euler formula on candidates |
| `field_assign.jl` | 152 | `_count_field_assignments_expanded()`, `_assign_expanded!()`, `_count_closed_loops_expanded()` | Backtracking field assignment with canonical multi-edge ordering; fermion-loop closure detection (union-find); ÷2^n for self-loop symmetry |
| `vertex_check.jl` | 73 | `_check_vertex_expanded()`, `_partial_vertex_ok_expanded()`, `_is_submultiset()` | Sorted-multiset vertex validation (expanded field conjugates); sub-multiset pruning for early failure detection |

---

## Cross-Cutting Observations

### 1. Hard-Coded vs Algorithmic Processes

**Hard-Coded (Stable):**
- **QED, EW, QCD model tables:** Field lists and vertex rules explicitly enumerated. No parser. Design choice: FeynRules anti-pattern avoided; rules are compile-time data, not runtime dicts.
- **Tree 2→2 channels:** Topology hardcoded to s/t/u. Filtering only: checks `rules.vertices[permutation]` for each channel. No qgraf call.

**Algorithmic (Phase 17c → Phase 18a-7):**
- **Diagram generation:** Nogueira's adjacency-matrix method (qgraf 4.0.6 §3).
  - `_degree_partitions()` → `_enumerate_topologies()` → `_fill_entry!()` with backtracking.
  - Dedup: `_is_canonical_topo()` delegates to `QgrafPort.is_canonical_feynman` (qg21 full lex-next, not pairwise swaps).
  - **Phase 17c (Session 24):** `count_diagrams()` now routes to `QgrafPort.count_diagrams_qg21` (Strategy C port). Legacy in-Julia version preserved as `_count_diagrams_legacy()` for regression testing.

### 2. Three-Model Relationship

| Model | Fermions | Bosons | Gauge | Use Case | Status |
|-------|----------|--------|-------|----------|--------|
| **QED** (1-3 gen) | e, μ, τ | γ | U(1) | e⁺e⁻→μ⁺μ⁻ (full pipeline: SPIRAL 0-7 complete) | Stable. Tree diagrams validated. |
| **EW** | e, νₑ | γ, Z, W | SU(2)×U(1) | e⁺e⁻→W⁺W⁻ (3 diagrams: s-γ, s-Z, t-νₑ) | Scaffolding. Full pipeline incomplete; Burnside orbit deferred Phase 18b. |
| **QCD** | q (color × 3) | g (8 colors) | SU(3) | gg→gg, qq̄→gg (test topology; color traces TBD) | Scaffolding. Color-trace contraction not yet in pipeline. |
| **φ³** | none | none | trivial | Test topology enumeration (pure graph combinatorics) | Test harness. All diagrams valid. |

### 3. Coupling Between Layer 3 and qgraf Port

**Dependency Chain (Layers 1-3 → qgraf):**
1. **Model** → `model_fields()`, `feynman_rules()` (rules dict).
2. **Layer 3 expansion:** `_expand_model_for_diagen()` splits particle/antiparticle (qgraf link-array style).
3. **Degree partition:** `_model_vertex_degrees(rules)` extracts allowed arities (3, 4, ...).
4. **Topology enum:** `_enumerate_topologies()` builds adjacency matrices.
5. **Canonical check:** `_is_canonical_topo()` → **`QgrafPort.is_canonical_feynman()`** (qg21 port, not in-Julia).
6. **Field assignment:** `_count_field_assignments_expanded()` validates vertices via `vertex_rules` from expanded model.

**Key Coupling Point:** Layer 3 algebraic methods (partition + topo enum + field assign) are **100% in-Julia**. The only qgraf port call from Layers 1-3 is:
- `QgrafPort.is_canonical_feynman()` — dedup check (lex-next permutation from qg21 §13156-13291).

**Data Flow (count_diagrams):**
```
Model + Rules 
  ↓ (expand: particle/antiparticle)
ExpandedModel (all_fields, conjugate, vertex_rules)
  ↓ (extract arities)
{degree set}
  ↓ (enumerate partitions)
{DegreePartition[]}
  ↓ (enumerate topologies × field assignments)
{FeynmanTopology + field bindings}
  → count_diagrams_qg21(strategy C)  [Phase 17c only]
```

### 4. Tech Debt & Deferrals

**Marked Deferrals (Search: "Phase 18"):**

| Issue | Location | Phase | Impact |
|-------|----------|-------|--------|
| Full algorithmic channel generation | `diagram_gen.jl:86` | TODO | Tree 2→2 still uses `tree_channels()` hard-coded; Phase 17c delegates to qgraf for general case. |
| 4-gluon vertex Lorentz structure | `rules.jl:93` | 18b-1 | Arity > 3 vertex factors not yet implemented. QCD 4-vertex exists but Lorentz is TBD. |
| Multi-vertex fermion lines | `qgraf/fermion_line.jl:65,77,81` | 18b | Compton tree, fermion loops need multi-step traversal. Phase 18a only handles single-line-per-external. |
| Boson polarisation | `qgraf/vertex_assemble.jl:114` | 18b (spin-sum) | Explicit polarisation tensor TBD; placeholder alg(1) used. |
| Multi-edge closed-loop Burnside | `qgraf/burnside_combine.jl:65-74` | 18b-3 | 3+ fermion-line bundles not yet supported. |

**No @test_broken annotations in Layers 1-3** (deferred work is in comments, not code).

### 5. Design Patterns

**Strengths:**
- **Holy traits** (Massive, Massless, Charged, Neutral) for orthogonal properties.
- **Type-parametric species dispatch** (Field{S} + vertex_structure on FieldSpecies).
- **Callable struct** FeynmanRules for natural `rules(field_names...)` interface.
- **Expanded model abstraction** (ExpandedModel) decouples fermion parity from rule table.

**Weaknesses (CLAUDE.md §24):**
- `MomentumSum` constructor returns `Union{Nothing, Momentum, MomentumSum}` (type instability).
- `gamma_pair` returns `Union{Nothing, Number, Pair}` (Algebra layer, not surveyed here).
- Global mutable state `_COLOUR_DUMMY_COUNTER` in colour layer (not surveyed).

### 6. Pipeline Coverage

Only **e⁺e⁻→μ⁺μ⁻** (QED tree-level) runs the full Layers 1-3-4-5-6 pipeline. All other processes:
- **e⁺e⁻→W⁺W⁻** (EW): Stops at Layer 3 (channels generated); Phase 18a-7 assembles amplitude but Burnside orbit (Phase 18b) incomplete.
- **Diagram enumeration (general):** Phase 17c delegates to qgraf port; Layer 3 in-Julia methods used for regression only.
- **QCD, φ³:** Topology only, no Lorentz structure or color trace.

### 7. Spiral Coverage

**Spirals 0-7** (Phases 1-13, Sessions 1-23): e⁺e⁻→μ⁺μ⁻ pipeline complete.  
**Spiral 8** (Phase 14, Session 24+): Bug fixes (type instabilities), γ₅ handling, Eps contraction.  
**Spiral 9** (Phase 17-18): Diagram generation + EW model integration.

---

## Summary

Layers 1-2 are **clean, well-tested, and stable**. Layer 3 is **under active development** with a two-track strategy:
1. **In-Julia methods** (partition/topo/field) for counting and regression.
2. **qgraf port** (Strategy C, Phase 17c) for production diagram generation.

The separation is clear: algebraic topology enumeration stays in Julia; canonical form checking delegates to the ported qg21 logic. Hard-coded processes (model tables, tree channels) are intentional scaffolding, not anti-patterns — they serve as reference implementations and will be gradually replaced by algorithmic generation as Phases 18a-18b complete.

**Key Gap:** Diagram generation is the strategic bottleneck. Once Phase 18a-7 (emission_to_amplitude) is fully tested, Spirals 9-10 can integrate EW and loop-level processes into the pipeline.
