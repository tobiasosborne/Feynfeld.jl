#  Phase 12b: vlis / vmap / lmap construction from xg (extension to qg10).
#
#  After qg21 emits a canonical topology (xg adjacency matrix), qg10 builds
#  three label arrays needed by qgen:
#
#    vlis[i]   = vertex visited at position i (canonical traversal order)
#    invlis[v] = inverse: position of vertex v in vlis
#    vmap[v,k] = the k-th neighbor of vertex v (in vlis order)
#    lmap[v,k] = back-pointer: vmap[vmap[v,k], lmap[v,k]] == v
#    rdeg[v]   = "real" degree at the time vertex v is visited
#    sdeg[v]   = rdeg[v] + self-loop count
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12028-12102.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, qg21_enumerate!,
                              compute_qg10_labels

@testset "Phase 12b: qg10 labels (vlis/vmap/lmap)" begin

    @testset "phi3 tree φ→φφ — single internal sink" begin
        # n=4, externals 1..3 connect to internal vertex 4.
        # Only one Step C topology emitted by qg21.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        snapshots = TopoState[]
        qg21_enumerate!(s) do state
            push!(snapshots, deepcopy(state))
        end
        @test length(snapshots) == 1
        labels = compute_qg10_labels(snapshots[1])

        @test labels.vlis[1:4] == Int8[1, 2, 3, 4]
        # Externals each have one neighbor: vertex 4.
        @test labels.vmap[1, 1] == Int8(4)
        @test labels.vmap[2, 1] == Int8(4)
        @test labels.vmap[3, 1] == Int8(4)
        # Internal vertex 4: three neighbors, in vlis order = (1, 2, 3).
        @test labels.vmap[4, 1] == Int8(1)
        @test labels.vmap[4, 2] == Int8(2)
        @test labels.vmap[4, 3] == Int8(3)
        # rdeg[4] = xn[4] = 3, sdeg = rdeg + gam[4,4] = 3 + 0
        @test labels.rdeg[4] == Int8(3)
        @test labels.sdeg[4] == Int8(3)
    end

    @testset "phi3 tree φφ→φφ — single internal edge" begin
        # n=6, vertex 5 connects to externals {1,2} + vertex 6;
        # vertex 6 connects to externals {3,4} + vertex 5.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        snapshots = TopoState[]
        qg21_enumerate!(s) do state
            push!(snapshots, deepcopy(state))
        end
        @test length(snapshots) == 1
        labels = compute_qg10_labels(snapshots[1])

        @test labels.vlis[1:6] == Int8[1, 2, 3, 4, 5, 6]
        @test labels.vmap[5, 1] == Int8(1)
        @test labels.vmap[5, 2] == Int8(2)
        @test labels.vmap[5, 3] == Int8(6)
        @test labels.vmap[6, 1] == Int8(3)
        @test labels.vmap[6, 2] == Int8(4)
        @test labels.vmap[6, 3] == Int8(5)
    end

    @testset "lmap is a self-consistent back-pointer (qg10:12113-12130)" begin
        # For every (v, k): vmap[vmap[v,k], lmap[v,k]] == v.
        # Vacuum diagrams (n_ext=0) excluded — qgraf errors on vaux=0
        # (see qg10:12055-12059 'qg10_1') and Feynfeld doesn't run them.
        for p in [
            Partition(Int8(3), Int8[1], Int8(3), Int8(0)),     # phi3 tree φ→φφ
            Partition(Int8(4), Int8[2], Int8(3), Int8(0)),     # phi3 tree φφ→φφ
        ]
            s = TopoState(p)
            qg21_enumerate!(s) do state
                labels = compute_qg10_labels(state)
                n = Int(state.n)
                for v in 1:n, k in 1:Int(state.vdeg[v])
                    j1 = Int(labels.vmap[v, k])
                    j2 = Int(labels.lmap[v, k])
                    @test labels.vmap[j1, j2] == Int8(v)
                end
            end
        end
    end

end
