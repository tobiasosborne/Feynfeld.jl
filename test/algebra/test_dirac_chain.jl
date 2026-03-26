# Tests for DiracChain (DOT) and DotSimplify
# Ref: FeynCalc NonCommAlgebra/DotSimplify.m

@testset "dot construction" begin
    # Basic chain
    chain = dot(GA(:μ), GA(:ν))
    @test chain isa DiracChain
    @test length(chain) == 2

    # Flattening
    c1 = dot(GA(:μ), GA(:ν))
    c2 = dot(GA(:ρ), GA(:σ))
    c3 = dot(c1, c2)
    @test length(c3) == 4

    # With spinors
    u = Spinor(:ubar, Momentum(:p), :m)
    v = Spinor(:u, Momentum(:k), :m)
    chain = dot(u, GA(:μ), v)
    @test length(chain) == 3
    @test chain.elements[1] isa Spinor
    @test chain.elements[3] isa Spinor
end

@testset "dot via * operator" begin
    chain = GA(:μ) * GA(:ν) * GA(:ρ)
    @test length(chain) == 3
    @test chain.elements[1].slot.index.name == :μ
    @test chain.elements[3].slot.index.name == :ρ
end

@testset "dot with slashed momenta" begin
    chain = dot(GA(:μ), GS(:p), GA(:ν))
    @test length(chain) == 3
    @test chain.elements[2].slot isa MomSlot
    @test chain.elements[2].slot.mom == Momentum(:p)
end

@testset "dot_simplify: no expansion needed" begin
    chain = dot(GA(:μ), GA(:ν))
    result = dot_simplify(chain)
    @test length(result) == 1
    @test result[1] == (1 // 1, chain)
end

@testset "dot_simplify: expand MomentumSum" begin
    # γ^μ (p+q)-slash → γ^μ p-slash + γ^μ q-slash
    pq = Momentum(:p) + Momentum(:q)
    chain = dot(GA(:μ), dirac_gamma(pq))
    result = dot_simplify(chain)
    @test length(result) == 2

    # Check both terms
    chains = Dict(length(r[2].elements) => r for r in result)
    for (_, (coeff, ch)) in chains
        @test coeff == 1 // 1
        @test length(ch) == 2
        @test ch.elements[1] == GA(:μ)
        @test ch.elements[2].slot isa MomSlot
        @test ch.elements[2].slot.mom isa Momentum
    end
end

@testset "dot_simplify: coefficients" begin
    # (2p - q)-slash → 2 * p-slash + (-1) * q-slash
    mom = 2 * Momentum(:p) - Momentum(:q)
    chain = DiracChain([dirac_gamma(mom)])
    result = dot_simplify(chain)
    @test length(result) == 2
    coeffs = Dict(r[2].elements[1].slot.mom.name => r[1] for r in result)
    @test coeffs[:p] == 2 // 1
    @test coeffs[:q] == -1 // 1
end

@testset "dot_simplify: multiple MomentumSums" begin
    # (p+q)-slash * (k+l)-slash → 4 terms
    pq = Momentum(:p) + Momentum(:q)
    kl = Momentum(:k) + Momentum(:l)
    chain = dot(dirac_gamma(pq), dirac_gamma(kl))
    result = dot_simplify(chain)
    @test length(result) == 4
end
