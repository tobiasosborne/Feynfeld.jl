#  Phase 16: momentum routing (spanning tree + chord identification).
#  Source: qg21:13315-13558.  Full leaf-peeling momentum assignment is a
#  Layer-4 concern handled when the new pipeline drives amplitude eval.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                              qg21_enumerate!, build_spanning_tree, count_chords

@testset "Phase 16: spanning tree + chord count" begin

    @testset "phi3 tree φ→φφ — 0 chords (tree topology)" begin
        # 4 vertices, 3 edges, no cycles ⇒ all edges in tree, 0 chords.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        @test count_chords(emits[1]) == 0
    end

    @testset "phi3 tree φφ→φφ — 0 chords (tree topology)" begin
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        @test count_chords(emits[1]) == 0
    end

    @testset "build_spanning_tree: phi3 tree φφ→φφ — 5 tree edges" begin
        # 6 vertices, 5 edges, all in tree.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        in_tree = build_spanning_tree(emits[1])
        @test count(in_tree) == 5
    end

    @testset "Euler invariant: chord count = nloop for non-vacuum" begin
        # For each emitted topology, chord count should equal the prescribed nloop.
        # Vacuum cases excluded (no externals → BFS from vertex 1 ok but
        # spanning tree visits all internals; chord count via Euler = nloop).
        for (n_ext, counts, mrho, nloop) in [
            (Int8(3), Int8[1], Int8(3), Int8(0)),                 # tree φ→φφ
            (Int8(4), Int8[2], Int8(3), Int8(0)),                 # tree φφ→φφ
        ]
            p = Partition(n_ext, counts, mrho, nloop)
            s = TopoState(p)
            qg21_enumerate!(s) do state
                @test count_chords(state) == Int(nloop)
            end
        end
    end

end
