#  Phase 17 prep: enumerate_topology_automorphisms.
#  Foundation for the dedup audition (Burnside / canonical-pmap / pre-filter).

using Test
using Feynfeld.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              enumerate_topology_automorphisms

@testset "enumerate_topology_automorphisms" begin

    @testset "phi3 tree φ→φφ — full auto group is 3! = 6 (ext perms)" begin
        # 1 internal vertex; 3 externals all attach to the same internal.
        # Swapping any 2 externals preserves gam → S_3 action on externals.
        # |G| = 3! = 6.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        autos = enumerate_topology_automorphisms(emits[1])
        @test length(autos) == 6
        @test autos[1] == Int8[1, 2, 3, 4]      # identity always first
    end

    @testset "phi3 tree φφ→φφ — 2 autos (identity + swap of 2 internals)" begin
        # 2 internals (5, 6) joined by 1 edge; symmetric.
        # Externals (1,2) attached to 5; (3,4) to 6 ⇒ external swap pairs
        # 1↔2, 3↔4 are NOT in the topology auto group (they don't change xg
        # because externals are vdeg=1; but they're in the equivalence class
        # of vdeg=1 vertices and our enumerator does iterate them).
        # Actually compute_equiv_classes! groups by (vdeg, xn, xg_diag) and
        # within that the perm enumerator runs them — so swaps of equivalent
        # externals appear iff they preserve gam.
        # Verify by counting:
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        autos = enumerate_topology_automorphisms(emits[1])
        # Topology: ext{1,2} connect to vertex 5; ext{3,4} to 6; v5↔v6.
        # Topology autos (preserving gam):
        #   identity
        #   swap (1,2)        — externals at vertex 5 swap, gam unchanged
        #   swap (3,4)        — externals at vertex 6 swap
        #   swap (1,2) ∧ (3,4)
        #   swap (5,6) ∧ (1↔3, 2↔4) — full mirror
        #   swap (5,6) ∧ (1↔4, 2↔3)
        #   swap (5,6) ∧ (1↔3, 2↔4) ∧ (1,2) — etc.
        # Total 8 autos under (vdeg, xn, xg_diag) classes. The actual count
        # depends on equiv-class grouping; the relation is:
        #   |G| = 2 (internal swap) × 2! × 2! (within-class ext swaps when
        #          pairs go to same internal) = 8.
        @test length(autos) == 8
        # Identity must be first.
        @test autos[1] == Int8[1, 2, 3, 4, 5, 6]
    end

    @testset "phi3 1L tadpole shape (1 ext) — 2 autos" begin
        # 1 ext + 2 deg-3 internals + 1L: half-edges = 1+6=7 → odd ✗
        # 3 ext + 3 deg-3 internals + 1L: 3+9=12 → P=6, V=6, L=1 ✓
        p = Partition(Int8(3), Int8[3], Int8(3), Int8(1))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        # Multiple topologies; just check enumerate runs and includes identity.
        for emit in emits
            autos = enumerate_topology_automorphisms(emit)
            @test length(autos) >= 1
            n = Int(emit.n)
            @test autos[1] == Int8.(1:n)
            # ngsym set as side-effect of compute_equiv_classes!; should equal length.
            @test Int(emit.ngsym) == length(autos)
        end
    end

end
