#  Phase 11: qg10 — external-leg permutation enumeration.
#  Source: qgraf-4.0.6.f08:12001-12200.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, qg10_enumerate!, _lex_next_perm!

@testset "qg10: external-leg permutation enumeration" begin

    @testset "_lex_next_perm!: Knuth Algorithm L" begin
        # All 6 permutations of [1,2,3] in lex order.
        p = Int8[1, 2, 3]
        seen = [copy(p)]
        while _lex_next_perm!(p)
            push!(seen, copy(p))
        end
        @test length(seen) == 6
        @test seen[1] == Int8[1, 2, 3]
        @test seen[2] == Int8[1, 3, 2]
        @test seen[3] == Int8[2, 1, 3]
        @test seen[4] == Int8[2, 3, 1]
        @test seen[5] == Int8[3, 1, 2]
        @test seen[6] == Int8[3, 2, 1]
    end

    @testset "_lex_next_perm!: singleton + empty" begin
        @test _lex_next_perm!(Int8[1]) == false
        @test _lex_next_perm!(Int8[]) == false
    end

    @testset "qg10: phi3 tree φ→φφ — 1 topo × 3! = 6 emissions" begin
        # 1 abstract topology × 6 ext-leg perms = 6 (qg10 raw, before qgen dedup).
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        emits = Tuple{Matrix{Int8}, Vector{Int8}}[]
        qg10_enumerate!(s) do state, perm
            n = Int(state.n)
            push!(emits, (copy(state.xg[1:n, 1:n]), copy(perm)))
        end
        @test length(emits) == 6
        # All ext-leg perms should be distinct.
        perms = unique([e[2] for e in emits])
        @test length(perms) == 6
    end

    @testset "qg10: phi3 tree φφ→φφ — 1 topo × 4! = 24 emissions" begin
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        nemit = Ref(0)
        qg10_enumerate!(s) do _, _
            nemit[] += 1
        end
        @test nemit[] == 24
        # Note: qgen (Phase 12) will dedupe these to 3 (s/t/u channels).
    end

    @testset "qg10: vacuum dumbbell — 2 topos × 0! = 2 emissions" begin
        # 0 ext, 2 int: ext perm is empty (just emit each topology once).
        p = Partition(Int8(0), Int8[2], Int8(3), Int8(2))
        s = TopoState(p)
        nemit = Ref(0)
        qg10_enumerate!(s) do _, perm
            nemit[] += 1
            @test isempty(perm)
        end
        @test nemit[] == 2
    end

end
