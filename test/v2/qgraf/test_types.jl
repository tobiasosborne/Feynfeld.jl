# Unit tests for qgraf port types.
# Beads: feynfeld-5hr

using Test
using Feynfeld.QgrafPort

@testset "qgraf port: types.jl" begin

    @testset "Partition: phi3 tree phi→phiphi" begin
        # 3 external legs, 1 internal 3-vertex, 0 loops
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        @test p.n_ext == 3
        @test p.nloop == 0
        @test rho_k(p, 3) == 1
        @test rho_k(p, 4) == 0
        @test rho_k(p, 2) == 0        # below mrho
        @test n_internal(p) == 1
        @test n_vertices(p) == 4      # 3 ext + 1 internal
        @test n_edges(p) == 3          # (3*1 + 3)/2 = 3 half-edges / 2
    end

    @testset "Partition: phi3 2-loop φφ → φφ (THE regression case)" begin
        # The Strategy C port target: 465 diagrams, qgraf says so.
        # Vertex accounting: Σ(k-2)ρ(k) = n_ext + 2(L-1)
        #   4 external legs, 2 loops  →  target = 4 + 2 = 6
        #   pure φ³: all internal vertices are degree 3, (3-2)·ρ(3) = 6  →  ρ(3)=6
        p = Partition(Int8(4), Int8[6], Int8(3), Int8(2))
        @test n_internal(p) == 6
        @test n_vertices(p) == 10
        # Σ k·ρ(k) + n_ext = 3·6 + 4 = 22  →  P = 11 edges
        @test n_edges(p) == 11
        # Euler: L = P − V + 1 = 11 − 10 + 1 = 2  ✓
        @test n_edges(p) - n_vertices(p) + 1 == p.nloop
    end

    @testset "Partition: mixed degrees (φ³ + φ⁴)" begin
        # Pretend-model with both 3- and 4-vertices.
        # counts[1] = ρ(3), counts[2] = ρ(4)
        p = Partition(Int8(2), Int8[2, 1], Int8(3), Int8(1))
        @test rho_k(p, 3) == 2
        @test rho_k(p, 4) == 1
        @test rho_k(p, 5) == 0
        @test n_internal(p) == 3
        @test n_vertices(p) == 5
        # Σ k·ρ(k) + n_ext = 3·2 + 4·1 + 2 = 12  →  P = 6
        @test n_edges(p) == 6
        # L = 6 - 5 + 1 = 2  (this partition has TWO loops, not 1)
        # Construct with nloop=2 to be consistent — this tests we don't
        # validate consistency at construction (validation happens in qpg11).
    end

    @testset "EquivClass: length + iteration" begin
        c = EquivClass(Int8(5), Int8(8), Int8(3))
        @test first(c) == 5
        @test last(c) == 8
        @test c.degree == 3
        @test length(c) == 4
        @test collect(c) == Int8[5, 6, 7, 8]

        singleton = EquivClass(Int8(3), Int8(3), Int8(4))
        @test length(singleton) == 1
        @test collect(singleton) == Int8[3]

        # Non-contiguous class: {3, 5, 7}
        sparse = EquivClass(Int8[3, 5, 7], Int8(3))
        @test length(sparse) == 3
        @test collect(sparse) == Int8[3, 5, 7]
        @test first(sparse) == 3
        @test last(sparse) == 7
    end

    @testset "FilterSet: default and keyword construction" begin
        # Default: everything off → no_filters is true
        f = FilterSet()
        @test no_filters(f)
        @test !f.onepi
        @test !f.nosnail
        @test !f.nosigma

        # One-hot construction
        f1 = FilterSet(; onepi=true)
        @test !no_filters(f1)
        @test f1.onepi
        @test !f1.nosnail

        # Multiple flags
        f2 = FilterSet(; onepi=true, nosnail=true, nosigma=true)
        @test f2.onepi
        @test f2.nosnail
        @test f2.nosigma
        @test !f2.onevi
    end

    @testset "TopoState: construction from Partition (phi3 tree)" begin
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        @test s.n == 4
        @test s.n_ext == 3
        @test s.rhop1 == 4          # first internal-vertex index
        @test s.nloop == 0

        # Vertex degrees: externals get 1, internal gets 3.
        @test s.vdeg[1] == 1
        @test s.vdeg[2] == 1
        @test s.vdeg[3] == 1
        @test s.vdeg[4] == 3

        # Scratch arrays sized to MAX_V, zero-initialised.
        @test length(s.vdeg) == MAX_V
        @test length(s.xn) == MAX_V
        @test size(s.xg) == (MAX_V, MAX_V)
        @test size(s.ds) == (MAX_V, MAX_V)
        @test all(iszero, s.xg)
        @test all(iszero, s.xn)
        @test all(iszero, s.xc)
        @test all(iszero, s.dta)
        @test s.ngsym == 0
        @test isempty(s.classes)
    end

    @testset "TopoState: phi3 2-loop shape" begin
        p = Partition(Int8(4), Int8[6], Int8(3), Int8(2))
        s = TopoState(p)
        @test s.n == 10
        @test s.n_ext == 4
        @test s.rhop1 == 5
        @test s.nloop == 2
        for i in 1:4
            @test s.vdeg[i] == 1
        end
        for i in 5:10
            @test s.vdeg[i] == 3
        end
        # Beyond n, vdeg stays zero (unused).
        for i in 11:MAX_V
            @test s.vdeg[i] == 0
        end
    end

    @testset "TopoState: mixed degree partition" begin
        # 2 external legs, ρ(3)=2, ρ(4)=1
        p = Partition(Int8(2), Int8[2, 1], Int8(3), Int8(2))
        s = TopoState(p)
        @test s.n == 5
        @test s.n_ext == 2
        @test s.rhop1 == 3
        # Internals ascend by degree: two 3-vertices then one 4-vertex.
        @test s.vdeg[3] == 3
        @test s.vdeg[4] == 3
        @test s.vdeg[5] == 4
    end

    @testset "TopoState: rejects n > MAX_V" begin
        # Construct a partition that would need more than MAX_V vertices.
        huge = Partition(Int8(0), Int8[MAX_V + 1], Int8(3), Int8(0))
        @test_throws Exception TopoState(huge)
    end

    @testset "EquivClass: rejects first > last + 1" begin
        # Empty class is allowed (first == last + 1 conventionally).
        @test_throws Exception EquivClass(Int8(5), Int8(3), Int8(3))
        # Singleton and empty are allowed:
        @test EquivClass(Int8(3), Int8(3), Int8(4)) isa EquivClass
        @test length(EquivClass(Int8(4), Int8(3), Int8(3))) == 0
    end

    @testset "FilterSet: all flags on" begin
        f = FilterSet(; onepi=true, nobridge=true, nosbridge=true, notadpole=true,
                       onshell=true, nosnail=true, onevi=true, onshellx=true,
                       noselfloop=true, nodiloop=true, noparallel=true,
                       cycli=true, nosigma=true, bipart=true)
        @test !no_filters(f)
        @test f.cycli && f.bipart && f.nosigma
    end
end
