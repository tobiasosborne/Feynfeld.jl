#  Phase 6: qg21 Step B — external-leg distribution (xc, xn enumeration).
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12554-12658
#
#  Tests are added incrementally.  Order:
#    Test 1: _degree_class_bounds helper (qg21:12541-12550)
#    Test 2: empty externals (rho(-1)==0) is a no-op (qg21:12554-12557)
#    Test 3: φ³ 1-loop tadpole emits exactly 1 (xc, xn) configuration
#    Test 4: QED ee→μμ tree emits exactly 5 configurations
#    Test 5: xn within each class is non-increasing (canonicalisation)

using Test
using Feynfeld.QgrafPort: Partition, TopoState, MAX_V, step_b_enumerate!
using Feynfeld.QgrafPort: _degree_class_bounds   # internal helper

@testset "qg21 Step B: external-leg distribution" begin

    @testset "_degree_class_bounds: phi3 tree φ→φφ" begin
        # vdeg = [1, 1, 1, 3]   (3 externals, 1 internal of degree 3)
        # Source: qg21:12541-12550 — xl[c]..xt[c] = contiguous run of equal vdeg.
        # Class 1 = {1,2,3} (vdeg=1); Class 2 = {4} (vdeg=3).
        vdeg = Int8[1, 1, 1, 3]
        xl, xt = _degree_class_bounds(vdeg, 4)
        @test xl == Int8[1, 4]
        @test xt == Int8[3, 4]
    end

    @testset "_degree_class_bounds: QED ee→μμ tree" begin
        # vdeg = [1, 1, 1, 1, 3, 3]   (4 externals, 2 internals of degree 3)
        vdeg = Int8[1, 1, 1, 1, 3, 3]
        xl, xt = _degree_class_bounds(vdeg, 6)
        @test xl == Int8[1, 5]
        @test xt == Int8[4, 6]
    end

    @testset "_degree_class_bounds: mixed φ³+φ⁴" begin
        # vdeg = [1, 1, 3, 3, 4]
        vdeg = Int8[1, 1, 3, 3, 4]
        xl, xt = _degree_class_bounds(vdeg, 5)
        @test xl == Int8[1, 3, 5]
        @test xt == Int8[2, 4, 5]
    end

    @testset "_degree_class_bounds: n=0 is empty" begin
        vdeg = Int8[]
        xl, xt = _degree_class_bounds(vdeg, 0)
        @test isempty(xl)
        @test isempty(xt)
    end

    # ─── enumeration core ─────────────────────────────────────────────────

    "Test helper: collect all (xc-prefix, xn-prefix) emissions from Step B."
    function _collect_step_b(state::TopoState)
        results = Tuple{Vector{Int8}, Vector{Int8}}[]
        n = Int(state.n)
        # xc has length n_classes+1 where n_classes is the number of internal
        # degree classes; cap at n+1 for simplicity.
        step_b_enumerate!(state) do s
            push!(results, (copy(s.xc[1:n+1]), copy(s.xn[1:n])))
        end
        results
    end

    @testset "Step B: empty externals — single no-op emission" begin
        # qg21:12554-12557:  if rho(-1)==0, set xc[1]=0 and emit once.
        # Vacuum 2-loop dumbbell: 2 deg-3 internals, 0 externals, nloop=2.
        p = Partition(Int8(0), Int8[2], Int8(3), Int8(2))
        s = TopoState(p)
        emits = _collect_step_b(s)
        @test length(emits) == 1
        @test s.xc[1] == Int8(0)
    end

    @testset "Step B: φ³ 1-loop tadpole — exactly 1 emission" begin
        # 1 ext, 2 int(deg-3), nloop=1.
        # Only one canonical xn distribution: xn[2]=1, xn[3]=0.
        # Hand-traced from qg21:12569-12658.
        p = Partition(Int8(1), Int8[2], Int8(3), Int8(1))
        s = TopoState(p)
        emits = _collect_step_b(s)
        @test length(emits) == 1
        # The single emission: xc=[0,1], xn=[0,1,0].
        xc, xn = emits[1]
        @test xc[1:2] == Int8[0, 1]
        @test xn       == Int8[0, 1, 0]
    end

    @testset "Step B: QED ee→μμ tree — exactly 5 emissions" begin
        # 4 ext, 2 int(deg-3), nloop=0.
        # Hand-traced enumeration:
        #   1. xc=[0,4], xn=[0,0,0,0,3,1]
        #   2. xc=[0,4], xn=[0,0,0,0,2,2]
        #   3. xc=[2,2], xn=[0,0,0,0,2,0]
        #   4. xc=[2,2], xn=[0,0,0,0,1,1]
        #   5. xc=[4,0], xn=[0,0,0,0,0,0]
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = _collect_step_b(s)
        @test length(emits) == 5
    end

    @testset "Step B: xn within each class is non-increasing" begin
        # Canonicalisation: xn[i] ≥ xn[i+1] for vertices in the same class.
        # qg21:12576-12583 — greedy non-increasing fill; backtrack preserves it.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = _collect_step_b(s)
        for (xc, xn) in emits
            # Internal class is vertices 5..6 (degree 3 each).
            @test xn[5] >= xn[6]
        end
    end

    @testset "Step B: multi-class partition (φ³+φ⁴ 2L, n_ext=2) — sanity" begin
        # 2 ext + 2 deg-3 + 1 deg-4 internal; 2 loops.
        # Half-edges = 2 + 6 + 4 = 12 → P=6, V=5, L=2.
        # vdeg = [1,1,3,3,4]; classes: ext (1..2), deg-3 (3..4), deg-4 (5).
        p = Partition(Int8(2), Int8[2, 1], Int8(3), Int8(2))
        s = TopoState(p)
        emits = _collect_step_b(s)

        # At least one emission and all emissions satisfy:
        #   (a) xc[1] is even and ≤ n_ext;
        #   (b) per-class xn sums equal xc per class;
        #   (c) per-class xn is non-increasing.
        @test length(emits) >= 1
        for (xc, xn) in emits
            @test xc[1] % Int8(2) == Int8(0)
            @test xc[1] <= Int8(2)
            # Class deg-3 (vertices 3..4): sum equals xc[2]; non-increasing.
            @test xn[3] + xn[4] == xc[2]
            @test xn[3] >= xn[4]
            # Class deg-4 (vertex 5): sum equals xc[3].
            @test xn[5] == xc[3]
        end
    end

end
