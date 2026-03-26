# Tests for Phase 1b remaining: DiracOrder, DiracEquation, DiracSimplify, Scheme
import Feynfeld: Pair

# ── DiracOrder ───────────────────────────────────────────────────────

@testset "DiracOrder: already ordered" begin
    chain = dot(GA(:μ), GA(:ν))  # μ < ν, already ordered
    result = dirac_order(chain)
    @test length(result) == 1
    @test result[1][1] == 1
    @test result[1][2] == chain
end

@testset "DiracOrder: swap two gammas" begin
    # g^ν g^μ → -g^μ g^ν + 2g^{μν}
    chain = dot(GA(:ν), GA(:μ))  # ν > μ, out of order
    result = dirac_order(chain)
    @test length(result) == 2
    # One term: -1 * ordered chain
    ordered_term = filter(r -> length(r[2]) == 2, result)
    @test length(ordered_term) == 1
    @test ordered_term[1][1] == -1
    # Other term: 2*metric * empty chain
    metric_term = filter(r -> length(r[2]) == 0, result)
    @test length(metric_term) == 1
end

@testset "DiracOrder: non-index gammas preserved" begin
    # Slashed momenta and g5 don't participate in ordering
    chain = dot(GA5(), GS(:p))
    result = dirac_order(chain)
    @test length(result) == 1
    @test result[1][2] == chain  # unchanged
end

# ── DiracEquation ────────────────────────────────────────────────────

@testset "DiracEquation: right boundary" begin
    # p-slash u(p) = m u(p)
    u = Spinor(:u, Momentum(:p), :m)
    chain = dot(GS(:p), u)
    result = dirac_equation(chain)
    @test length(result) == 1
    @test result[1][1] == :m  # coefficient is mass
    @test result[1][2].elements == [u]
end

@testset "DiracEquation: right boundary v-spinor" begin
    # p-slash v(p) = -m v(p)
    v = Spinor(:v, Momentum(:p), :m)
    chain = dot(GS(:p), v)
    result = dirac_equation(chain)
    @test length(result) == 1
    c = result[1][1]
    # _mul_coeff(-1, :m) produces :(-1 * m) or just -:m
    @test c isa Expr || c == -1  # symbolic -m
end

@testset "DiracEquation: left boundary" begin
    # ubar(p) p-slash = m ubar(p)
    ubar = Spinor(:ubar, Momentum(:p), :m)
    chain = dot(ubar, GS(:p))
    result = dirac_equation(chain)
    @test length(result) == 1
    @test result[1][1] == :m
    @test result[1][2].elements == [ubar]
end

@testset "DiracEquation: no match" begin
    # g^μ u(p) — gamma is not p-slash, no simplification
    u = Spinor(:u, Momentum(:p), :m)
    chain = dot(GA(:μ), u)
    result = dirac_equation(chain)
    @test result[1][1] == 1
    @test result[1][2] == chain
end

# ── DiracSimplify ────────────────────────────────────────────────────

@testset "DiracSimplify: g^mu g_mu" begin
    result = dirac_simplify(dot(GA(:μ), GA(:μ)))
    @test length(result) == 1
    @test result[1][1] == 4
    @test length(result[1][2]) == 0
end

@testset "DiracSimplify: p-slash squared" begin
    result = dirac_simplify(dot(GS(:p), GS(:p)))
    @test length(result) == 1
    @test result[1][1] == SP(:p, :p)
end

@testset "DiracSimplify: expand + simplify" begin
    # (p+q)-slash . (p+q)-slash expands to 4 terms:
    # p.p-slash + p.q-slash + q.p-slash + q.q-slash
    # p-slash.p-slash → SP(:p,:p), q-slash.q-slash → SP(:q,:q)
    # p-slash.q-slash and q-slash.p-slash → chains of 2 (different momenta)
    pq = Momentum(:p) + Momentum(:q)
    chain = dot(dirac_gamma(pq), dirac_gamma(pq))
    result = dirac_simplify(chain)
    nonzero = filter(r -> r[1] != 0, result)
    @test length(nonzero) >= 2
    # The p^2 and q^2 terms have empty chains
    scalar_terms = filter(r -> length(r[2]) == 0, nonzero)
    @test length(scalar_terms) >= 2
end

@testset "DiracSimplify: with Dirac equation" begin
    u = Spinor(:u, Momentum(:p), :m)
    ubar = Spinor(:ubar, Momentum(:p), :m)
    # ubar(p) p-slash u(p) = m * ubar(p) u(p) → m^2 * chain
    # Actually: ubar . p-slash . u(p) → left: m*ubar . u(p)
    chain = dot(ubar, GS(:p), u)
    result = dirac_simplify(chain)
    @test length(result) >= 1
    # Should have simplified the p-slash away
    @test any(r -> length(r[2]) == 2, result)  # ubar, u remaining
end

# ── DiracScheme ──────────────────────────────────────────────────────

@testset "Scheme system" begin
    @test current_scheme() == NDR

    set_scheme!(BMHV_SCHEME)
    @test current_scheme() == BMHV_SCHEME

    # with_scheme restores
    with_scheme(LARIN) do
        @test current_scheme() == LARIN
    end
    @test current_scheme() == BMHV_SCHEME

    # Restore default
    set_scheme!(NDR)
    @test current_scheme() == NDR
end
