# Phase 1b validation: translated FeynCalc Dirac MUnit tests
# Source: refs/FeynCalc/Tests/Dirac/DiracTrace.test, DiracSimplify.test
import Feynfeld: Pair

# ── DiracTrace MUnit translations ────────────────────────────────────

@testset "MUnit Tr[1] = 4" begin
    @test dirac_trace(DiracChain(DiracElement[])) == 4
end

@testset "MUnit Tr[g^a g^b] = 4 g^{ab}" begin
    # fcstDiracTrace-ID: Tr[GA[a].GA[b]] = 4 MT[a,b]
    result = dirac_trace(dot(GA(:a), GA(:b)))
    @test result == Feynfeld._mul_coeff(4, MT(:a, :b))
end

@testset "MUnit Tr[odd gammas] = 0" begin
    @test dirac_trace(dot(GA(:a))) == 0
    @test dirac_trace(dot(GA(:a), GA(:b), GA(:c))) == 0
    @test dirac_trace(dot(GS(:p))) == 0
end

@testset "MUnit Tr[g5] = 0" begin
    @test dirac_trace(dot(GA5())) == 0
end

@testset "MUnit Tr[g^a g^b g^c g^d] (4-gamma trace)" begin
    # fcstDiracTrace: Tr[GA[a].GA[b].GA[c].GA[d]]
    # = 4(g^{ab}g^{cd} - g^{ac}g^{bd} + g^{ad}g^{bc})
    result = dirac_trace(dot(GA(:a), GA(:b), GA(:c), GA(:d)))
    @test result isa Vector
    @test length(result) == 3  # three terms
end

@testset "MUnit Tr[p-slash q-slash] = 4 p.q" begin
    result = dirac_trace(dot(GS(:p), GS(:q)))
    @test result == Feynfeld._mul_coeff(4, SP(:p, :q))
end

@testset "MUnit Tr[p-slash q-slash k-slash l-slash]" begin
    # Tr[GS[p].GS[q].GS[k].GS[l]] = 4(p.q k.l - p.k q.l + p.l q.k)
    result = dirac_trace(dot(GS(:p), GS(:q), GS(:k), GS(:l)))
    @test result isa Vector
    @test length(result) == 3
end

@testset "MUnit Tr[g5 g^a g^b g^c g^d] (chiral trace)" begin
    # Tr[GA5().GA[a].GA[b].GA[c].GA[d]] = -4 * Eps[a,b,c,d]
    result = dirac_trace(dot(GA5(), GA(:a), GA(:b), GA(:c), GA(:d)))
    @test result != 0
end

@testset "MUnit Tr[g5 g^a g^b] = 0" begin
    @test dirac_trace(dot(GA5(), GA(:a), GA(:b))) == 0
end

# ── DiracTrick MUnit translations ────────────────────────────────────

@testset "MUnit g^mu g_mu = 4 (4D)" begin
    result = dirac_trick(dot(GA(:mu), GA(:mu)))
    @test result[1][1] == 4
end

@testset "MUnit g^mu g_mu = D (D-dim)" begin
    result = dirac_trick(dot(GAD(:mu), GAD(:mu)))
    @test result[1][1] == :D
end

@testset "MUnit p-slash p-slash = p^2" begin
    result = dirac_trick(dot(GS(:p), GS(:p)))
    @test result[1][1] == SP(:p, :p)
end

@testset "MUnit g^mu g^nu g_mu = -2 g^nu (4D)" begin
    result = dirac_trick(dot(GA(:mu), GA(:nu), GA(:mu)))
    @test result[1][1] == -2
    @test result[1][2].elements == [GA(:nu)]
end

@testset "MUnit g5 g5 = 1" begin
    result = dirac_trick(dot(GA5(), GA5()))
    @test result[1][1] == 1
    @test length(result[1][2]) == 0
end

@testset "MUnit g6 g7 = 0" begin
    result = dirac_trick(dot(GA6(), GA7()))
    @test all(r -> r[1] == 0, result)
end

@testset "MUnit g6 g6 = g6" begin
    result = dirac_trick(dot(GA6(), GA6()))
    @test result[1][2].elements == [GA6()]
end

# ── DiracSimplify MUnit translations ─────────────────────────────────

@testset "MUnit DiracSimplify g^mu g_mu" begin
    result = dirac_simplify(dot(GA(:mu), GA(:mu)))
    @test result[1][1] == 4
end

@testset "MUnit DiracSimplify p-slash p-slash" begin
    result = dirac_simplify(dot(GS(:p), GS(:p)))
    @test result[1][1] == SP(:p, :p)
end

# ── DiracOrder MUnit translations ────────────────────────────────────

@testset "MUnit DiracOrder: g^b g^a → -g^a g^b + 2 g^{ab}" begin
    result = dirac_order(dot(GA(:b), GA(:a)))
    @test length(result) == 2
    # One term with ordered chain, one with metric
    ordered = filter(r -> length(r[2]) == 2, result)
    @test length(ordered) == 1
    @test ordered[1][1] == -1  # sign from anticommutation
    metric = filter(r -> length(r[2]) == 0, result)
    @test length(metric) == 1
end

# ── DiracEquation MUnit translations ─────────────────────────────────

@testset "MUnit DiracEquation: p-slash u(p) = m u(p)" begin
    u = Spinor(:u, Momentum(:p), :m)
    chain = dot(GS(:p), u)
    result = dirac_equation(chain)
    @test result[1][1] == :m
    @test result[1][2].elements == [u]
end
