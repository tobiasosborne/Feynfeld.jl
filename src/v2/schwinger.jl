# Schwinger correction: O(α) QED correction to σ(e+e- → μ+μ-).
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 7.2 + Section 8
# Cross-check: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/El-GaEl.m
# "δσ/σ = (3α)/(4π) × (π²/3 - 1/2)"
#
# The total NLO correction (virtual + soft real, inclusive) in the
# high-energy limit s >> m_e², m_μ² is:
#   δσ/σ = (3α)/(4π) × (π²/3 - 1/2)
#
# This is the vertex + soft real emission correction only (VP not included).
# IR divergences cancel between virtual and real emission (KLN theorem).

# ──────────────────────────────────────────────────────────────────
# REFERENCE IMPLEMENTATION — validates physics formula directly.
# The pipeline should reproduce this result via
# Model → Rules → Diagrams → Algebra → Integrals → Evaluate.
# Do NOT delete: used for cross-validation against pipeline results.
# ──────────────────────────────────────────────────────────────────

using QuadGK: quadgk

"""
    schwinger_correction(; alpha=1/137.036)

O(α) QED correction to σ(e⁺e⁻ → μ⁺μ⁻) total cross-section.
Returns δσ/σ in the high-energy limit (s >> m²).
"""
schwinger_correction(; alpha::Float64 = 1 / 137.036) =
    3 * alpha / (4 * π) * (π^2 / 3 - 1 / 2)

"""
    vacuum_polarization(s, m2; alpha, mu2)

Renormalized vacuum polarization Π̂(s) = Π(s) - Π(0) for one fermion flavour
with mass² = m2. On-shell scheme: Π̂(0) = 0.

Computed via Feynman parameter integral using QuadGK:
  Π̂(s) = -(2α/π) ∫₀¹ dx x(1-x) ln[1 - x(1-x)s/(m² - iε)]
"""
function vacuum_polarization(s::Float64, m2::Float64;
                              alpha::Float64 = 1 / 137.036)::ComplexF64
    s == 0.0 && return 0.0 + 0.0im
    val, _ = quadgk(0.0, 1.0; rtol = 1e-12) do x
        arg = complex(1.0 - x * (1.0 - x) * s / m2, 1e-30)
        x * (1.0 - x) * log(arg)
    end
    -2.0 * alpha / π * val
end

"""
    sigma_nlo_ee_mumu(s; alpha)

NLO total cross-section σ(e⁺e⁻ → μ⁺μ⁻) including O(α) Schwinger correction.
"""
function sigma_nlo_ee_mumu(s::Float64; alpha::Float64 = 1 / 137.036)
    sigma_0 = 4 * π * alpha^2 / (3 * s)
    sigma_0 * (1.0 + schwinger_correction(; alpha))
end
