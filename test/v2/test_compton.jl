# Spiral 1: Compton scattering e(p1) + γ(k1) → e(p2) + γ(k2)
#
# Ground truth: Peskin & Schroeder Eq. (5.87)
# Cited from: refs/FeynCalc/FeynCalc/Examples/QED/Tree/Mathematica/ElGa-ElGa.m, lines 103-108
# "|M̄|² = 2e⁴ [ SP[p1,k2]/SP[p1,k1] + SP[p1,k1]/SP[p1,k2]
#              + 2m² (1/SP[p1,k1] - 1/SP[p1,k2])
#              + m⁴ (1/SP[p1,k1] - 1/SP[p1,k2])² ]"
#
# Two diagrams: s-channel + u-channel (both added, same sign).
# Ref: refs/FeynCalc/.../ElGa-ElGa.md confirms "CORRECT" vs P&S 5.87.

using Test
include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

@testset "Spiral 1: Compton scattering" begin

    @testset "PolarizationSum" begin
        # Source: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID1
        # "PolarizationSum[rho, si] = -Pair[LorentzIndex[rho], LorentzIndex[si]]"
        mu = LorentzIndex(:mu)
        nu = LorentzIndex(:nu)
        @test polarization_sum(mu, nu) == -alg(pair(mu, nu))
    end

    @testset "Compton |M̄|² from algebra pipeline" begin
        # Momenta
        p1 = Momentum(:p1)
        p2 = Momentum(:p2)
        k1 = Momentum(:k1)
        k2 = Momentum(:k2)

        # The two Compton diagrams contribute gamma chains between spinors:
        # Γ₁ = γ^ν (p̸₁+k̸₁+m) γ^μ   with denominator d₁ = s - m²
        # Γ₂ = γ^μ (p̸₁-k̸₂+m) γ^ν   with denominator d₂ = u - m²
        #
        # After spin sum + polarisation sum:
        # |M̄|² = (e⁴/4) × Σ_{i,j} (-g_{μμ'})(-g_{νν'}) Tr[(p̸₂+m)Γᵢ(p̸₁+m)Γ̃ⱼ] / (dᵢ dⱼ)
        #
        # Γ̃ⱼ = reverse of Γⱼ with indices relabeled: μ→μ', ν→ν'
        # The (-g)(-g) from pol sums contracts μ-μ' and ν-ν'.

        # Build each trace numerically at specific kinematics.
        # Choose: m² = 1, s = 5, u = -2, t = 2m²-s-u = -1
        m2_val = 1//1
        s_val = 5//1
        u_val = -2//1
        t_val = 2 * m2_val - s_val - u_val  # = -1

        # Scalar products from Mandelstam:
        # p1² = m², p2² = m², k1² = 0, k2² = 0
        # p1·k1 = (s-m²)/2, p1·k2 = (m²-u)/2
        # p2·k2 = (s-m²)/2, p2·k1 = (m²-u)/2  (crossing symmetry)
        # p1·p2 = (2m²-t)/2 = (2+1)/2 = 3/2
        # k1·k2 = -t/2 = 1/2
        p1k1 = (s_val - m2_val) // 2  # = 2
        p1k2 = (m2_val - u_val) // 2  # = 3/2

        ctx = sp_context(
            (:p1, :p1) => m2_val,
            (:p2, :p2) => m2_val,
            (:k1, :k1) => 0//1,
            (:k2, :k2) => 0//1,
            (:p1, :k1) => p1k1,
            (:p1, :k2) => p1k2,
            (:p2, :k1) => p1k2,     # crossing
            (:p2, :k2) => p1k1,     # crossing
            (:p1, :p2) => (2 * m2_val - t_val) // 2,
            (:k1, :k2) => -t_val // 2,
        )

        # Propagator denominators
        d1 = s_val - m2_val  # = 4
        d2 = u_val - m2_val  # = -3

        # For each (i,j) diagram pair, build the full trace gamma chain:
        # Tr[(p̸₂+m) Γᵢ (p̸₁+m) Γ̃ⱼ]
        #
        # We expand (p̸+m) as p-slash + m·identity, then sum over the 4 combinations.
        # The trace of a chain with identity γ is just the trace of the remaining gammas.
        #
        # Helper: compute trace for a single (i,j) pair
        function compton_trace_ij(gamma_i, gamma_j_conj, p_out, m_out, p_in, m_in)
            # Tr[(p̸_out + m_out) × gamma_i × (p̸_in + m_in) × gamma_j_conj]
            # = Tr[p̸_out × Γᵢ × p̸_in × Γ̃ⱼ]
            # + m_in × Tr[p̸_out × Γᵢ × Γ̃ⱼ]
            # + m_out × Tr[Γᵢ × p̸_in × Γ̃ⱼ]
            # + m_in × m_out × Tr[Γᵢ × Γ̃ⱼ]
            result = AlgSum()

            # Term 1: Tr[p̸_out Γᵢ p̸_in Γ̃ⱼ]
            gs1 = [GS(p_out); gamma_i; GS(p_in); gamma_j_conj]
            result = result + dirac_trace(gs1)

            # Term 2: m_in × Tr[p̸_out Γᵢ Γ̃ⱼ]
            if !iszero(m_in)
                gs2 = [GS(p_out); gamma_i; gamma_j_conj]
                result = result + m_in * dirac_trace(gs2)
            end

            # Term 3: m_out × Tr[Γᵢ p̸_in Γ̃ⱼ]
            if !iszero(m_out)
                gs3 = [gamma_i; GS(p_in); gamma_j_conj]
                result = result + m_out * dirac_trace(gs3)
            end

            # Term 4: m_in × m_out × Tr[Γᵢ Γ̃ⱼ]
            if !iszero(m_in) && !iszero(m_out)
                gs4 = [gamma_i; gamma_j_conj]
                result = result + m_in * m_out * dirac_trace(gs4)
            end

            result
        end

        # Gamma chains for the two diagrams (between completeness insertions):
        # Γ₁ = γ^ν (p̸₁+k̸₁) γ^μ + m γ^ν γ^μ
        #     → propagator numerator p̸₁+k̸₁+m as sum of chains
        # BUT: we need to separate the propagator numerator from m.
        # The trace handles MomentumSum expansion.

        # Build p1+k1 as MomentumSum
        p1pk1 = MomentumSum([(1//1, p1), (1//1, k1)])
        p1mk2 = MomentumSum([(1//1, p1), (-1//1, k2)])

        # Γ₁ gammas (propagator numerator included, mass separate):
        # Full: γ^ν (p̸₁+k̸₁+m) γ^μ
        # = γ^ν (p̸₁+k̸₁) γ^μ + m γ^ν γ^μ
        gamma_1_mom = DiracGamma[GAD(:nu), DiracGamma(MomSumSlot(p1pk1)), GAD(:mu)]
        gamma_1_mass = DiracGamma[GAD(:nu), GAD(:mu)]

        # Γ₂ gammas:
        # Full: γ^μ (p̸₁-k̸₂+m) γ^ν
        gamma_2_mom = DiracGamma[GAD(:mu), DiracGamma(MomSumSlot(p1mk2)), GAD(:nu)]
        gamma_2_mass = DiracGamma[GAD(:mu), GAD(:nu)]

        # Conjugate gammas (reversed order, relabeled indices):
        # Γ̃₁: γ^{μ'} (p̸₁+k̸₁) γ^{ν'} + m γ^{μ'} γ^{ν'}
        gamma_1c_mom = DiracGamma[GAD(:mu_), DiracGamma(MomSumSlot(p1pk1)), GAD(:nu_)]
        gamma_1c_mass = DiracGamma[GAD(:mu_), GAD(:nu_)]

        # Γ̃₂: γ^{ν'} (p̸₁-k̸₂) γ^{μ'} + m γ^{ν'} γ^{μ'}
        gamma_2c_mom = DiracGamma[GAD(:nu_), DiracGamma(MomSumSlot(p1mk2)), GAD(:mu_)]
        gamma_2c_mass = DiracGamma[GAD(:nu_), GAD(:mu_)]

        m_val = m2_val  # mass = 1 (m² = 1, so m = 1 as rational)

        # For each (i,j), compute Tr[(p̸₂+m)Γᵢ(p̸₁+m)Γ̃ⱼ] by expanding
        # Γᵢ = Γᵢ_mom + m × Γᵢ_mass and similar for Γ̃ⱼ.
        # This gives 4 sub-traces per (i,j) from the propagator numerator expansion.

        function trace_diagram_pair(gi_mom, gi_mass, gjc_mom, gjc_mass)
            result = AlgSum()
            # momentum × momentum
            result = result + compton_trace_ij(gi_mom, gjc_mom, p2, m_val, p1, m_val)
            # momentum × mass (conjugate)
            result = result + m_val * compton_trace_ij(gi_mom, gjc_mass, p2, m_val, p1, m_val)
            # mass × momentum (conjugate)
            result = result + m_val * compton_trace_ij(gi_mass, gjc_mom, p2, m_val, p1, m_val)
            # mass × mass (conjugate)
            result = result + m_val^2 * compton_trace_ij(gi_mass, gjc_mass, p2, m_val, p1, m_val)
            result
        end

        # Compute T₁₁, T₁₂, T₂₁, T₂₂
        T11 = trace_diagram_pair(gamma_1_mom, gamma_1_mass, gamma_1c_mom, gamma_1c_mass)
        T12 = trace_diagram_pair(gamma_1_mom, gamma_1_mass, gamma_2c_mom, gamma_2c_mass)
        T21 = trace_diagram_pair(gamma_2_mom, gamma_2_mass, gamma_1c_mom, gamma_1c_mass)
        T22 = trace_diagram_pair(gamma_2_mom, gamma_2_mass, gamma_2c_mom, gamma_2c_mass)

        # Contract Lorentz indices (μ-μ' and ν-ν' from pol sums, plus internal)
        # The pol sums give factors of (-g^{μμ'})(-g^{νν'}) = g^{μμ'} g^{νν'}
        # When we contract, μ-μ' and ν-ν' get contracted.
        # But we also have the metric from the photon propagator in the trace.

        # Actually: the pol sum replaces ε_μ ε*_{μ'} → -g_{μμ'}.
        # In the squared amplitude, both photon pol sums give (-g)(-g).
        # The two (-1) factors give +1 overall for the polarisation contribution.
        # The contraction of μ with μ' and ν with ν' is done by `contract`.

        # Combine with denominators and spin average (1/4 for 2 initial spins × 2 initial pols):
        # |M̄|² = (1/4) × [T₁₁/(d₁²) + T₁₂/(d₁ d₂) + T₂₁/(d₂ d₁) + T₂₂/(d₂²)]
        # Ref: refs/FeynCalc/.../ElGa-ElGa.m line 88: ExtraFactor -> 1/2^2
        total = 1//4 * ((1 // d1^2) * T11 + (1 // (d1 * d2)) * T12 +
                        (1 // (d2 * d1)) * T21 + (1 // d2^2) * T22)

        # Apply polarization sums: relabel primed → unprimed indices.
        # Feynman gauge: Σ_pol ε^μ ε^{μ'*} = -g^{μμ'} → contract μ with μ'.
        # Relabeling μ'→μ and ν'→ν makes them repeated; contract handles the rest.
        # Factor: (-1)² = +1 from two polarization sums (already included).
        # NB: Must use DimD() to match GAD-created indices.
        mu_  = LorentzIndex(:mu_, DimD())
        nu_  = LorentzIndex(:nu_, DimD())
        mu_idx = LorentzIndex(:mu, DimD())
        nu_idx = LorentzIndex(:nu, DimD())
        pol_summed = substitute_index(total, mu_, mu_idx)
        pol_summed = substitute_index(pol_summed, nu_, nu_idx)

        # Contract all Lorentz indices
        contracted = contract(pol_summed; ctx)

        # Expand scalar products (MomentumSum bilinear expansion)
        expanded = expand_scalar_product(contracted)

        # Evaluate all scalar products at our kinematics
        result = evaluate_sp(expanded; ctx)

        # Extract the numerical value (should be a pure number)
        @test length(result.terms) == 1  # single scalar term
        fk, coeff = first(result.terms)
        @test isempty(fk.factors)  # no remaining factors

        pipeline_value = coeff isa DimPoly ? evaluate_dim(coeff) : coeff

        # Compare to P&S Eq. 5.87
        # Ref: refs/FeynCalc/.../ElGa-ElGa.m, lines 103-105
        # "|M̄|²/(2e⁴) = SP[p1,k2]/SP[p1,k1] + SP[p1,k1]/SP[p1,k2]
        #              + 2m² (1/SP[p1,k1] - 1/SP[p1,k2])
        #              + m⁴ (1/SP[p1,k1] - 1/SP[p1,k2])²"
        ps587 = 2 * (p1k2 // p1k1 + p1k1 // p1k2 +
                     2 * m2_val * (1 // p1k1 - 1 // p1k2) +
                     m2_val^2 * (1 // p1k1 - 1 // p1k2)^2)

        # The (1/4) spin average is already included in our calculation.
        # The e⁴ coupling is stripped (we compute the reduced matrix element).
        @test pipeline_value == ps587
    end
end
