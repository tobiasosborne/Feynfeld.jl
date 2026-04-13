#  Phase 3: 4-gluon vertex in qcd_model.
#
#  qgraf's QCD model includes a quartic gluon vertex.
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/models/qcd verbatim:
#     [ gluon, gluon, gluon, gluon ]
#  Lorentz structure: V^{abcd}_{μνρσ} ∝ f^{abe}f^{cde}(g_{μρ}g_{νσ} - g_{μσ}g_{νρ})
#                                       + 2 perms (Peskin & Schroeder Eq. 16.10).
#
#  This test only checks DIAGRAM COUNTS — the Lorentz tensor is a Layer-4
#  concern handled in a later phase.  The 4-gluon contact diagram closes
#  the previously-missing 1 + 1 + 10 diagrams across these three processes.

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "QCD 4-gluon vertex (golden master)" begin

    qcd = qcd_model()

    # gg → gg tree: 4 diagrams (s + t + u via ggg + 1 contact via gggg)
    # Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/golden_masters/qcd/
    #         gluon_gluon_TO_gluon_gluon_0L.out
    @test count_diagrams(qcd, [:g, :g], [:g, :g]; loops=0) == 4

    # qq̄ → ggg tree: 16 diagrams (15 with only ggg + 1 with gggg contact)
    # Source: golden_masters/qcd/quark_antiquark_TO_gluon_gluon_gluon_0L.out
    @test count_diagrams(qcd, [:q, :q], [:g, :g, :g]; loops=0) == 16

    # gg → ggg tree: 25 diagrams (15 with only ggg + 10 contact channels)
    # Source: golden_masters/qcd/gluon_gluon_TO_gluon_gluon_gluon_0L.out
    @test count_diagrams(qcd, [:g, :g], [:g, :g, :g]; loops=0) == 25

end
