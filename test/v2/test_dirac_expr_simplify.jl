# Regression: simplify(::DiracExpr) must NOT drop spinors when collecting
# like terms. Two terms whose chains differ ONLY in their spinor pair are
# physically distinct and must remain so after simplify (and therefore
# after `+`, which calls simplify).
#
# Bug fixed in feynfeld-ocpb (F002): pre-fix simplify keyed the term
# Dict on `gammas(chain)` (gamma-only filter), so two terms with the
# same gammas but different spinors collided and merged.

using Test
using Feynfeld

@testset "DiracExpr simplify preserves spinors" begin
    p1 = Momentum(:p1)
    p2 = Momentum(:p2)

    # Two chains with identical gammas (γ^μ) but different spinor pairs:
    #   chain_a = [ūbar(p1), γ^μ, v(p2)]   — line 1: out-electron + out-positron
    #   chain_b = [v̄(p2), γ^μ, u(p1)]      — line 2: in-positron + in-electron
    # These represent physically distinct fermion lines (e+e- vs annihilation
    # in opposite direction); they MUST NOT merge under simplify.
    chain_a = dot(ubar(p1), GA(:mu), v(p2))
    chain_b = dot(vbar(p2), GA(:mu), u(p1))

    de_a = DiracExpr([(alg(2), chain_a)])
    de_b = DiracExpr([(alg(3), chain_b)])

    # Direct construction: two terms preserved trivially
    de_direct = DiracExpr([(alg(2), chain_a), (alg(3), chain_b)])
    @test length(de_direct.terms) == 2

    # Through `+` (which calls simplify): terms must NOT collapse
    sum_ab = de_a + de_b
    @test length(sum_ab.terms) == 2

    # The full chains (with spinors) must survive
    chains_after = Set([t[2] for t in sum_ab.terms])
    @test chain_a in chains_after
    @test chain_b in chains_after

    # Coefficients must NOT have merged: each chain keeps its own coefficient
    coeff_a = first(t[1] for t in sum_ab.terms if t[2] == chain_a)
    coeff_b = first(t[1] for t in sum_ab.terms if t[2] == chain_b)
    @test coeff_a == alg(2)
    @test coeff_b == alg(3)
end

@testset "DiracExpr simplify: identical chains DO merge" begin
    p1 = Momentum(:p1)
    p2 = Momentum(:p2)

    # Two terms with identical chains (same gammas AND same spinors)
    # SHOULD merge their coefficients (this is what simplify is for).
    chain = dot(ubar(p1), GA(:mu), v(p2))
    de1 = DiracExpr([(alg(2), chain)])
    de2 = DiracExpr([(alg(3), chain)])
    sum = de1 + de2

    @test length(sum.terms) == 1
    @test sum.terms[1][1] == alg(5)   # 2 + 3
    @test sum.terms[1][2] == chain
end

@testset "DiracExpr simplify: spinor-free chains still merge correctly" begin
    # Chains with no spinors (typical of dirac_trick output): same gammas
    # SHOULD merge — this is the existing dirac_trick.jl pathway.
    chain_a = DiracChain(DiracElement[GA(:mu), GA(:nu)])
    chain_b = DiracChain(DiracElement[GA(:mu), GA(:nu)])
    de_a = DiracExpr([(alg(2), chain_a)])
    de_b = DiracExpr([(alg(3), chain_b)])
    sum = de_a + de_b
    @test length(sum.terms) == 1
    @test sum.terms[1][1] == alg(5)
end
