# Unit tests for qgraf canonical.jl — the 474/465 bug locus.
# Beads: feynfeld-xjc

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort
using .FeynfeldX.QgrafPort: _lex_next!, _compare_permuted_adjacency

# ─── Helper: build a bare TopoState for unit tests ────────────────────────
function _bare_state(n_ext, n, vdeg_internal)
    p = Partition(Int8(n_ext), Int8[length(vdeg_internal)], Int8(3), Int8(0))
    # p.counts is a placeholder; we override vdeg below.
    s = TopoState(p)
    s.n = Int8(n)
    s.n_ext = Int8(n_ext)
    s.rhop1 = Int8(n_ext + 1)
    fill!(s.vdeg, Int8(0))
    for i in 1:n_ext
        s.vdeg[i] = Int8(1)
    end
    for (k, d) in enumerate(vdeg_internal)
        s.vdeg[n_ext + k] = Int8(d)
    end
    s
end

@testset "qgraf port: canonical.jl" begin

    @testset "_lex_next!: Knuth Algorithm L" begin
        perm = Int8[1, 2, 3, 0, 0, 0, 0, 0]
        @test _lex_next!(perm, 1, 3) && perm[1:3] == Int8[1, 3, 2]
        @test _lex_next!(perm, 1, 3) && perm[1:3] == Int8[2, 1, 3]
        @test _lex_next!(perm, 1, 3) && perm[1:3] == Int8[2, 3, 1]
        @test _lex_next!(perm, 1, 3) && perm[1:3] == Int8[3, 1, 2]
        @test _lex_next!(perm, 1, 3) && perm[1:3] == Int8[3, 2, 1]
        @test !_lex_next!(perm, 1, 3)   # exhausted: 3,2,1 is max
    end

    @testset "_lex_next!: singleton and adjacent pair" begin
        perm = Int8[5, 5, 0, 0]
        @test !_lex_next!(perm, 1, 1)   # size 1: nothing to permute
        perm = Int8[1, 2, 0, 0]
        @test _lex_next!(perm, 1, 2) && perm[1:2] == Int8[2, 1]
        @test !_lex_next!(perm, 1, 2)
    end

    @testset "next_class_perm!: single class of size 3 = 3! permutations" begin
        perm = Int8[1, 2, 3, 0, 0, 0, 0, 0]
        for i in 1:3; perm[i] = Int8(i); end
        classes = [EquivClass(Int8(1), Int8(3), Int8(3))]
        seen = Set{NTuple{3, Int8}}()
        push!(seen, (perm[1], perm[2], perm[3]))
        while next_class_perm!(perm, classes)
            push!(seen, (perm[1], perm[2], perm[3]))
        end
        @test length(seen) == 6   # 3!
    end

    @testset "next_class_perm!: two classes = product enumeration" begin
        # class 1 = {3,4}, class 2 = {5,6,7}
        perm = zeros(Int8, 8)
        for i in 1:8; perm[i] = Int8(i); end
        classes = [
            EquivClass(Int8(3), Int8(4), Int8(3)),
            EquivClass(Int8(5), Int8(7), Int8(3)),
        ]
        seen = Set{NTuple{5, Int8}}()
        push!(seen, (perm[3], perm[4], perm[5], perm[6], perm[7]))
        while next_class_perm!(perm, classes)
            push!(seen, (perm[3], perm[4], perm[5], perm[6], perm[7]))
        end
        @test length(seen) == 2 * 6   # 2! × 3! = 12
    end

    @testset "compute_equiv_classes!: one uniform class" begin
        # 2 externals + 3 internals all (deg=3, xn=0, xg_diag=0) → one class
        s = _bare_state(2, 5, [3, 3, 3])
        @test compute_equiv_classes!(s) == 1
        c = s.classes[1]
        @test c.members == Int8[3, 4, 5]
        @test c.degree == 3
    end

    @testset "compute_equiv_classes!: split by xn" begin
        # vertex 3 has xn=2, vertex 4,5 have xn=0 → two classes
        s = _bare_state(2, 5, [3, 3, 3])
        s.xn[3] = Int8(2); s.xn[4] = Int8(0); s.xn[5] = Int8(0)
        n = compute_equiv_classes!(s)
        @test n == 2
        # one class {3}, one class {4,5} (order is deterministic by key)
        all_members = sort(vcat([collect(c.members) for c in s.classes]...))
        @test all_members == Int8[3, 4, 5]
        @test Int8[3] in [c.members for c in s.classes]
        @test Int8[4, 5] in [c.members for c in s.classes]
    end

    @testset "compute_equiv_classes!: split by xg_diag" begin
        # vertex 4 has a self-loop (xg[4,4]=2) → own class
        s = _bare_state(2, 5, [3, 3, 3])
        s.xg[4, 4] = Int8(2)
        n = compute_equiv_classes!(s)
        @test n == 2  # {3,5} share xg_diag=0, {4} alone
        all_members = sort(vcat([collect(c.members) for c in s.classes]...))
        @test all_members == Int8[3, 4, 5]
        @test Int8[3, 5] in [c.members for c in s.classes]
        @test Int8[4] in [c.members for c in s.classes]
    end

    @testset "compute_equiv_classes!: non-contiguous grouping" begin
        # vertices 3 and 5 share invariant; 4 is alone.
        # xn[3]=0, xn[4]=2, xn[5]=0 → {3,5} and {4}.
        s = _bare_state(2, 5, [3, 3, 3])
        s.xn[3] = Int8(0); s.xn[4] = Int8(2); s.xn[5] = Int8(0)
        n = compute_equiv_classes!(s)
        @test n == 2
        # The class {3,5} is non-contiguous — this was the bug case.
        class35 = [c for c in s.classes if c.members == Int8[3, 5]]
        class4 = [c for c in s.classes if c.members == Int8[4]]
        @test length(class35) == 1
        @test length(class4) == 1
    end

    @testset "compute_equiv_classes!: no internals" begin
        # only externals
        s = _bare_state(2, 2, Int[])
        @test compute_equiv_classes!(s) == 0
        @test isempty(s.classes)
    end

    @testset "is_canonical_full!: trivial (no internals)" begin
        s = _bare_state(2, 2, Int[])
        @test is_canonical_full!(s)
        @test s.ngsym == 1
    end

    @testset "is_canonical_full!: 2 equivalent internals joined by 1 edge → ngsym=2" begin
        # 2 ext + 2 internal, single edge 3-4, both in same class
        s = _bare_state(2, 4, [3, 3])
        s.xg[3, 4] = Int8(1)   # upper-triangular storage
        @test is_canonical_full!(s)
        @test s.ngsym == 2   # identity + swap
    end

    @testset "is_canonical_full!: swap yields smaller matrix → NOT canonical" begin
        # xg[3,4]=1 only, all else 0. All three internals {3,4,5} share
        # (vdeg=3, xn=0, xg_diag=0) → one class. Under perm swap (4,5) the
        # "1" moves from position (3,4) [earlier in row-major] to (3,5)
        # [later], so the permuted matrix is lex-SMALLER than the original
        # at the earlier position (3,4). Under lex-smallest convention,
        # the original is NOT canonical.
        s = _bare_state(2, 5, [3, 3, 3])
        s.xg[3, 4] = Int8(1)
        @test !is_canonical_full!(s)
    end

    @testset "is_canonical_full!: triangle (all 3 internals pairwise connected) → ngsym=6" begin
        # 0 ext + 3 internal, full triangle: each pair has one edge
        s = _bare_state(0, 3, [2, 2, 2])
        s.xg[1, 2] = Int8(1)
        s.xg[1, 3] = Int8(1)
        s.xg[2, 3] = Int8(1)
        s.rhop1 = Int8(1)   # no externals
        @test is_canonical_full!(s)
        @test s.ngsym == 6   # full S_3
    end

    @testset "_compare_permuted_adjacency: identity returns 0" begin
        s = _bare_state(2, 4, [3, 3])
        s.xg[3, 4] = Int8(1)
        compute_equiv_classes!(s)
        perm = Int8[1, 2, 3, 4, 0, 0, 0, 0]
        @test _compare_permuted_adjacency(s, perm) == 0
    end
end
