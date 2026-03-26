# Tests for AlgTerm/AlgSum expression tree
# Source: Feynfeld.jl/src/algebra/alg_expr.jl, alg_ops.jl

using Feynfeld: AlgTerm, AlgSum, AlgFactor, alg, alg_zero, alg_scalar, is_scalar,
    Pair, Eps, LorentzIndex, Momentum, MomentumSum, Dim4, DimD,
    SP, FV, MT, SPD, FVD, MTD, LC,
    DiracChain, DiracGamma, LISlot, MomSlot,
    dot, GA, GS,
    contract, expand_scalar_product, dirac_trace_alg,
    SPContext, set_sp, levi_civita

# ── AlgTerm construction ──────────────────────────────────────────

@testset "AlgTerm construction" begin
    t = AlgTerm(3)
    @test t.coeff == 3
    @test isempty(t.factors)

    p = SP(:p, :q)
    t2 = AlgTerm(1, p)
    @test t2.coeff == 1
    @test length(t2.factors) == 1
    @test t2.factors[1] == p

    t3 = AlgTerm(2, AlgFactor[MT(:μ, :ν), FV(:p, :ρ)])
    @test t3.coeff == 2
    @test length(t3.factors) == 2
end

# ── AlgSum construction and lifters ───────────────────────────────

@testset "AlgSum lifters" begin
    @test alg(42) isa AlgSum
    @test length(alg(42).terms) == 1
    @test alg(42).terms[1].coeff == 42

    @test alg(:D) isa AlgSum
    @test alg(:D).terms[1].coeff == :D

    p = SP(:p, :q)
    s = alg(p)
    @test s isa AlgSum
    @test length(s.terms) == 1
    @test s.terms[1].coeff == 1
    @test s.terms[1].factors[1] == p

    @test is_scalar(alg(5))
    @test !is_scalar(alg(MT(:μ, :ν)))

    @test isempty(alg_zero().terms)
end

# ── Arithmetic: addition ──────────────────────────────────────────

@testset "AlgSum addition" begin
    a = alg(SP(:p, :q))
    b = alg(SP(:k, :l))
    s = a + b
    @test s isa AlgSum
    @test length(s.terms) == 2

    # Adding zero
    @test length((a + alg_zero()).terms) == 1

    # Adding scalars
    s2 = alg(3) + alg(4)
    @test length(s2.terms) == 1
    @test s2.terms[1].coeff == 7

    # Collecting like terms: a + a = 2a
    doubled = a + a
    @test length(doubled.terms) == 1
    @test doubled.terms[1].coeff == 2
end

@testset "AlgSum subtraction" begin
    a = alg(SP(:p, :q))
    z = a - a
    @test isempty(z.terms)

    b = alg(3) - alg(1)
    @test length(b.terms) == 1
    @test b.terms[1].coeff == 2
end

# ── Arithmetic: multiplication ────────────────────────────────────

@testset "AlgSum multiplication" begin
    a = alg(MT(:μ, :ν))
    b = alg(FV(:p, :ρ))
    prod = a * b
    @test length(prod.terms) == 1
    @test prod.terms[1].coeff == 1
    @test length(prod.terms[1].factors) == 2

    # Scalar multiply
    s = 3 * alg(SP(:p, :q))
    @test s.terms[1].coeff == 3

    # FOIL: (a+b)(c+d) = ac + ad + bc + bd
    p1 = alg(MT(:μ, :ν)) + alg(FV(:k, :μ))
    p2 = alg(FV(:p, :ν)) + alg(SP(:p, :q))
    result = p1 * p2
    @test length(result.terms) == 4
end

@testset "AlgSum scalar-Pair multiply" begin
    s = alg(3) * alg(MT(:μ, :ν))
    @test length(s.terms) == 1
    @test s.terms[1].coeff == 3
    @test length(s.terms[1].factors) == 1
end

# ── contract(::AlgSum) ────────────────────────────────────────────

@testset "contract AlgSum" begin
    # g^{μν} g_{μρ} = g^{ν}_{ρ} ... expressed as AlgSum
    s = alg(MT(:μ, :ν)) * alg(MT(:μ, :ρ))
    c = contract(s)
    @test length(c.terms) == 1
    @test c.terms[1].coeff == 1
    @test c.terms[1].factors[1] == MT(:ν, :ρ)

    # g^{μ}_{μ} = 4
    s2 = alg(MT(:μ, :μ))
    c2 = contract(s2)
    @test is_scalar(c2)
    @test c2.terms[1].coeff == 4

    # p^μ q_μ = p·q
    s3 = alg(FV(:p, :μ)) * alg(FV(:q, :μ))
    c3 = contract(s3)
    @test length(c3.terms) == 1
    @test c3.terms[1].factors[1] == SP(:p, :q)

    # (g^{μν} + g^{μρ}) * p_ν = p^μ + contract(g^{μρ} p_ν) [no shared index]
    s4 = (alg(MT(:μ, :ν)) + alg(MT(:μ, :ρ))) * alg(FV(:p, :ν))
    c4 = contract(s4)
    # First term contracts to FV(:p, :μ), second has no shared index
    @test length(c4.terms) == 2
end

@testset "contract AlgSum with SPContext" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :q, :s)
    s = alg(FV(:p, :μ)) * alg(FV(:q, :μ))
    c = contract(s; ctx=ctx)
    @test is_scalar(c)
    @test c.terms[1].coeff == :s
end

# ── expand_scalar_product(::AlgSum) ───────────────────────────────

@testset "expand_scalar_product AlgSum" begin
    pq = Momentum(:p) + Momentum(:q)
    sp_pq_k = Feynfeld.Pair(pq, Momentum(:k))
    s = alg(sp_pq_k)
    expanded = expand_scalar_product(s)
    @test length(expanded.terms) == 2
    # Should contain SP(:k,:p) and SP(:k,:q)
    pairs_found = Set([t.factors[1] for t in expanded.terms])
    @test SP(:k, :p) in pairs_found
    @test SP(:k, :q) in pairs_found
end

@testset "expand_scalar_product AlgSum with context" begin
    ctx = SPContext()
    ctx = set_sp(ctx, :p, :p, :m2)
    pq = Momentum(:p) + Momentum(:q)
    sp = Feynfeld.Pair(pq, Momentum(:p))
    s = alg(sp)
    expanded = expand_scalar_product(s; ctx=ctx)
    # (p+q)·p = p² + q·p = m2 + SP(:p,:q)
    @test length(expanded.terms) == 2
end

# ── dirac_trace_alg ───────────────────────────────────────────────

@testset "dirac_trace_alg" begin
    # Tr[1] = 4
    empty_chain = DiracChain(Union{DiracGamma,Feynfeld.Spinor}[])
    r0 = dirac_trace_alg(empty_chain)
    @test is_scalar(r0)
    @test r0.terms[1].coeff == 4

    # Tr[g^μ g^ν] = 4 g^{μν}
    chain2 = dot(GA(:μ), GA(:ν))
    result = dirac_trace_alg(chain2)
    @test length(result.terms) == 1
    t = result.terms[1]
    @test t.coeff == 4
    @test length(t.factors) == 1
    @test t.factors[1] == MT(:μ, :ν)

    # Tr[g^μ g^ν g^ρ g^σ] = 4(g^{μν}g^{ρσ} - g^{μρ}g^{νσ} + g^{μσ}g^{νρ})
    chain4 = dot(GA(:μ), GA(:ν), GA(:ρ), GA(:σ))
    result4 = dirac_trace_alg(chain4)
    @test length(result4.terms) == 3
    # Each term should have coeff ±4 and exactly 2 Pair factors
    for t in result4.terms
        @test abs(t.coeff) == 4
        @test length(t.factors) == 2
    end
end

# ── Integration: trace × trace → contract ─────────────────────────

@testset "trace product contraction" begin
    # Core use case: Tr[g^μ g^ν] × Tr[g_μ g_ν]
    # = (4 g^{μν}) × (4 g^{μν}) = 16 g^{μν} g^{μν} = 16 * 4 = 64
    chain_a = dot(GA(:μ), GA(:ν))
    chain_b = dot(GA(:μ), GA(:ν))
    tr_a = dirac_trace_alg(chain_a)
    tr_b = dirac_trace_alg(chain_b)
    product = tr_a * tr_b
    contracted = contract(product)
    @test is_scalar(contracted)
    @test contracted.terms[1].coeff == 64
end

@testset "trace product with momenta" begin
    # Tr[p̸ g^μ] × Tr[q̸ g_μ] = (4 p^μ) × (4 q_μ) = 16 p·q
    chain_a = dot(GS(:p), GA(:μ))
    chain_b = dot(GS(:q), GA(:μ))
    tr_a = dirac_trace_alg(chain_a)
    tr_b = dirac_trace_alg(chain_b)
    product = tr_a * tr_b
    contracted = contract(product)
    @test length(contracted.terms) == 1
    t = contracted.terms[1]
    @test t.coeff == 16
    @test t.factors[1] == SP(:p, :q)
end
