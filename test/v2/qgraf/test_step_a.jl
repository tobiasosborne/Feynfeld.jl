#  Phase 5: qg21 Step A audit — degree-sequence init in TopoState.
#
#  qgraf source (Fortran 2008):
#    refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12479-12492
#  Verbatim:
#      n=0
#      do i1=1,rho(-1)
#        n=n+1
#        vdeg(n)=1
#      end do
#      jj=rho(-1)
#      do i1=mrho,nrho
#        jj=jj+i1*rho(i1)
#        do i2=1,rho(i1)
#          n=n+1
#          vdeg(n)=i1
#        end do
#      end do
#      nli=jj/2
#      loop=nli-n+1
#
#  TopoState(::Partition) reproduces this layout exactly.  These tests pin
#  the source-line correspondence and assert the Euler invariant
#  L = P − V + 1 on every partition shape that downstream phases will use.

using Test
using Feynfeld.QgrafPort: Partition, TopoState, n_vertices, n_internal,
                              n_edges, rho_k

@testset "qg21 Step A: degree-seq init (audit)" begin

    @testset "vdeg layout: phi3 tree φ → φφ" begin
        # n_ext=3 (1 in + 2 out), 1 internal vertex of degree 3, nloop=0.
        # Source line 12479: externals first with vdeg=1.
        # Source line 12485: internals next in ascending degree order.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        @test s.n        == Int8(4)
        @test s.n_ext    == Int8(3)
        @test s.rhop1    == Int8(4)
        @test s.vdeg[1:4] == Int8[1, 1, 1, 3]
    end

    @testset "vdeg layout: mixed φ³ + φ⁴" begin
        # n_ext=2, 1 internal of degree 3, 1 of degree 4. Internals ascend.
        p = Partition(Int8(2), Int8[1, 1], Int8(3), Int8(2))
        s = TopoState(p)
        @test s.n        == Int8(4)
        @test s.vdeg[1:4] == Int8[1, 1, 3, 4]
    end

    @testset "Euler invariant L = P − V + 1" begin
        # Source lines 12492-93: nli = jj/2; loop = nli − n + 1.
        # We reproduce this via Partition.nloop and n_edges/n_vertices.
        for (n_ext, counts, mrho, nloop) in [
            (Int8(3), Int8[1], Int8(3), Int8(0)),                 # tree φ→φφ
            (Int8(2), Int8[2], Int8(3), Int8(1)),                 # φφ→φφ 1L
            (Int8(2), Int8[4], Int8(3), Int8(2)),                 # φφ→φφ 2L
            (Int8(2), Int8[2, 1], Int8(3), Int8(2)),              # mixed: 2×deg-3 + 1×deg-4
                                                                    # half-edges = 2 + 3·2 + 4·1 = 12, P=6, V=5, L=2
            (Int8(4), Int8[2], Int8(3), Int8(0)),                 # φφ→φφ tree
        ]
            p = Partition(n_ext, counts, mrho, nloop)
            P = n_edges(p)
            V = n_vertices(p)
            L = P - V + 1
            @test L == Int(nloop)
        end
    end

    @testset "rejects n > MAX_V" begin
        p = Partition(Int8(20), Int8[10], Int8(3), Int8(0))   # n = 30
        @test_throws ErrorException TopoState(p)
    end

end
