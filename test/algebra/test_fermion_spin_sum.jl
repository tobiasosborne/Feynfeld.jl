# Tests for FermionSpinSum
# Source: Feynfeld.jl/src/algebra/fermion_spin_sum.jl
#
# Physics: completeness relations Σ_s u ū = p̸+m, Σ_s v v̄ = p̸-m
# Ref: Peskin & Schroeder §5.1

using Feynfeld: fermion_spin_sum, AlgSum, AlgTerm, alg, alg_zero, is_scalar,
    Spinor, DiracChain, DiracGamma, dot, GA, GS, GA5,
    Pair, SP, FV, MT, LorentzIndex, Momentum,
    contract, expand_scalar_product, dirac_trace_alg,
    SPContext, set_sp

# ── Basic completeness: massless ──────────────────────────────────

@testset "FermionSpinSum massless" begin
    # Simplest case: Σ [ū(p) γ^μ u(p)] [ū(p) γ^ν u(p)]
    # = Tr[p̸ γ^μ p̸ γ^ν] = 4(2 p^μ p^ν - p² g^{μν})
    chain_M = dot(Spinor(:ubar, Momentum(:p)), GA(:μ), Spinor(:u, Momentum(:p)))
    chain_Mc = dot(Spinor(:ubar, Momentum(:p)), GA(:ν), Spinor(:u, Momentum(:p)))
    result = fermion_spin_sum(chain_M, chain_Mc)
    @test result isa AlgSum
    @test !isempty(result.terms)
end

@testset "FermionSpinSum u-v pair massless" begin
    # e+e- → mu+mu- muon line (massless): Σ [ū(k) γ^μ v(k')] [v̄(k') γ^ν u(k)]
    # = Tr[k̸ γ^μ k̸' γ^ν]
    # = 4(k^μ k'^ν + k^ν k'^μ - k·k' g^{μν})
    k = Momentum(:k)
    kp = Momentum(:kp)
    chain_M = dot(Spinor(:ubar, k), GA(:μ), Spinor(:v, kp))
    chain_Mc = dot(Spinor(:vbar, kp), GA(:ν), Spinor(:u, k))
    result = fermion_spin_sum(chain_M, chain_Mc)
    @test result isa AlgSum
    # Should have 3 terms: 4*FV(k,μ)*FV(kp,ν) + 4*FV(k,ν)*FV(kp,μ) - 4*SP(k,kp)*MT(μ,ν)
    @test length(result.terms) == 3
    for t in result.terms
        @test abs(t.coeff) == 4
        @test length(t.factors) == 2
    end
end

# ── Massive case ──────────────────────────────────────────────────

@testset "FermionSpinSum massive u-u" begin
    # Σ [ū(p) u(p)] [ū(p) u(p)] with mass m
    # = Tr[(p̸+m)(p̸+m)] = Tr[p̸² + 2m p̸ + m²]
    # = p² Tr[1] + 2m Tr[p̸] + m² Tr[1]
    # = 4 p² + 0 + 4 m²  (Tr[p̸] = 0 for odd # gammas)
    p = Momentum(:p)
    chain_M = dot(Spinor(:ubar, p, :m), Spinor(:u, p, :m))
    chain_Mc = dot(Spinor(:ubar, p, :m), Spinor(:u, p, :m))
    result = fermion_spin_sum(chain_M, chain_Mc)
    @test result isa AlgSum
    # Should have terms: 4*SP(p,p) (from p̸ p̸ trace) and some m² terms
    @test !isempty(result.terms)
end

@testset "FermionSpinSum massive u-v pair" begin
    # Σ [ū(k) γ^μ v(k')] [v̄(k') γ^ν u(k)] with masses m_k, m_kp
    # = Tr[(k̸+m_k) γ^μ (k̸'-m_kp) γ^ν]
    # 4 terms from distributing mass factors
    k = Momentum(:k)
    kp = Momentum(:kp)
    chain_M = dot(Spinor(:ubar, k, :mk), GA(:μ), Spinor(:v, kp, :mkp))
    chain_Mc = dot(Spinor(:vbar, kp, :mkp), GA(:ν), Spinor(:u, k, :mk))
    result = fermion_spin_sum(chain_M, chain_Mc)
    @test result isa AlgSum
    # Massless part: 3 terms (from Tr[k̸ γ^μ k̸' γ^ν])
    # Mass corrections: additional terms from m_k and m_kp cross terms
    @test length(result.terms) >= 3
end

# ── Integration with contract ─────────────────────────────────────

@testset "FermionSpinSum + contract: two fermion lines" begin
    # Core e+e- → mu+mu- calculation (massless limit):
    # Line 1: Tr[k̸ γ^μ k̸' γ^ν]
    # Line 2: Tr[p̸ γ_ν p̸' γ_μ]
    # Product contracted over μ,ν
    k = Momentum(:k)
    kp = Momentum(:kp)
    p = Momentum(:p)
    pp = Momentum(:pp)

    # Muon line
    mu_M = dot(Spinor(:ubar, k), GA(:μ), Spinor(:v, kp))
    mu_Mc = dot(Spinor(:vbar, kp), GA(:ν), Spinor(:u, k))
    tr_mu = fermion_spin_sum(mu_M, mu_Mc)

    # Electron line
    e_M = dot(Spinor(:vbar, pp), GA(:μ), Spinor(:u, p))
    e_Mc = dot(Spinor(:ubar, p), GA(:ν), Spinor(:v, pp))
    tr_e = fermion_spin_sum(e_M, e_Mc)

    # Multiply and contract
    product = tr_mu * tr_e
    contracted = contract(product)
    @test contracted isa AlgSum
    @test !isempty(contracted.terms)
    # Should reduce to scalar products only (all Lorentz indices contracted)
    for t in contracted.terms
        for f in t.factors
            if f isa Pair
                @test !(f.a isa LorentzIndex && f.b isa LorentzIndex)
            end
        end
    end
end

# ── Error handling ────────────────────────────────────────────────

@testset "FermionSpinSum validation" begin
    k = Momentum(:k)
    kp = Momentum(:kp)
    # Momentum mismatch
    chain_M = dot(Spinor(:ubar, k), GA(:μ), Spinor(:v, kp))
    chain_bad = dot(Spinor(:vbar, k), GA(:ν), Spinor(:u, kp))  # wrong: vbar(k) should be vbar(kp)
    @test_throws ErrorException fermion_spin_sum(chain_M, chain_bad)
end
