# Tests for DiracTrick (core Dirac algebra simplification)
# Ref: FeynCalc Tests/Dirac/DiracTrick.test
import Feynfeld: Pair

@testset "gamma5 squared" begin
    # g^5 . g^5 = 1
    result = dirac_trick(dot(GA5(), GA5()))
    @test length(result) == 1
    @test result[1][1] == 1
    @test length(result[1][2]) == 0

    # g^5 . g^5 . g^mu = g^mu
    result = dirac_trick(dot(GA5(), GA5(), GA(:μ)))
    @test length(result) == 1
    @test result[1][1] == 1
    @test length(result[1][2]) == 1
    @test result[1][2].elements[1] == GA(:μ)
end

@testset "projector algebra" begin
    # g^6 . g^6 = g^6 (idempotent)
    result = dirac_trick(dot(GA6(), GA6()))
    @test length(result) == 1
    @test result[1][1] == 1
    @test result[1][2].elements == [GA6()]

    # g^7 . g^7 = g^7
    result = dirac_trick(dot(GA7(), GA7()))
    @test length(result) == 1
    @test result[1][2].elements == [GA7()]

    # g^6 . g^7 = 0
    result = dirac_trick(dot(GA6(), GA7()))
    @test all(r[1] == 0 for r in result)

    # g^7 . g^6 = 0
    result = dirac_trick(dot(GA7(), GA6()))
    @test all(r[1] == 0 for r in result)

    # g^5 . g^6 = g^6
    result = dirac_trick(dot(GA5(), GA6()))
    @test result[1][2].elements == [GA6()]

    # g^5 . g^7 = -g^7
    result = dirac_trick(dot(GA5(), GA7()))
    @test result[1][1] == -1
    @test result[1][2].elements == [GA7()]
end

@testset "slash squared" begin
    # p-slash . p-slash = p^2
    result = dirac_trick(dot(GS(:p), GS(:p)))
    @test length(result) == 1
    @test result[1][1] == SP(:p, :p)
    @test length(result[1][2]) == 0

    # Different momenta: no simplification
    result = dirac_trick(dot(GS(:p), GS(:q)))
    @test length(result) == 1
    @test length(result[1][2]) == 2
end

@testset "adjacent trace: g^mu g_mu = D" begin
    # 4D: g^mu g_mu = 4
    result = dirac_trick(dot(GA(:μ), GA(:μ)))
    @test result[1][1] == 4
    @test length(result[1][2]) == 0

    # D-dim: g^mu g_mu = D
    result = dirac_trick(dot(GAD(:μ), GAD(:μ)))
    @test result[1][1] == :D
    @test length(result[1][2]) == 0
end

@testset "sandwich: g^mu g^a g_mu = (2-D) g^a" begin
    # 4D: g^mu g^nu g_mu = -2 g^nu
    result = dirac_trick(dot(GA(:μ), GA(:ν), GA(:μ)))
    @test length(result) == 1
    @test result[1][1] == -2  # 2 - 4 = -2
    @test result[1][2].elements == [GA(:ν)]

    # D-dim: g^mu g^nu g_mu = (2-D) g^nu
    result = dirac_trick(dot(GAD(:μ), GAD(:ν), GAD(:μ)))
    @test length(result) == 1
    @test result[1][1] == :(2 - D)
    @test result[1][2].elements == [GAD(:ν)]
end

@testset "sandwich: g^mu g^a g^b g_mu" begin
    # 4D: g^mu g^a g^b g_mu = 4*g^{a,b} + (4-4)*g^a g^b = 4*g^{a,b}
    result = dirac_trick(dot(GA(:μ), GA(:α), GA(:β), GA(:μ)))
    nonzero = filter(r -> r[1] != 0, result)
    @test length(nonzero) == 1
    c, ch = nonzero[1]
    @test length(ch) == 0  # scalar result (empty chain)
    # Coefficient is 4*MT(:α,:β) — a product expression
    @test c isa Expr  # :(4 * Pair(...))
end
