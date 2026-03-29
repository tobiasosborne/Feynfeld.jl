# Running fine-structure constant α(q²) from vacuum polarization.
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Section 3, Eq. (3.10)
#      "α(s) = α / (1 - (Δα(s))_ferm)" [on-shell running coupling]
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eq. (10.11)-(10.12)
#      "α⁻¹(0) = 137.035999178(8)", "Δα_had^(5)(M_Z) = 0.02783 ± 0.00006"
#
# The renormalized vacuum polarization for fermion f (charge Q_f, N_c colors):
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eq. (4.23) + on-shell subtraction
# "Π̂(q²) = Π(q²) - Π(0)", with Π from B₀ via Eq. (4.23)
# In Feynman parameter form (standard textbook result, e.g. P&S Ch. 7):
#   Π̂_f(q²) = Q_f² N_c × [-(2α/π) ∫₀¹ dx x(1-x) ln(1 - x(1-x)q²/m_f²)]
#
# Running coupling (1-loop resummation):
#   α(q²) = α / (1 - Δα(q²)),   Δα = -Σ_f Π̂_f(q²)
#
# For timelike q² (s > 0), Π̂ is complex above threshold 4m_f². The physical
# cross-section uses |α(s)|² = α² / |1 - Π̂(s)|².
#
# Note: Light quark masses (u, d, s) enter perturbatively. Non-perturbative
# QCD effects make the hadronic contribution uncertain at ~10%. For precision
# work, use dispersion-relation data for Δα_had instead.

# ---- Standard Model fermion table (PDG 2024 values) ----
# Each entry: (mass² in GeV², charge², color factor N_c)
# Ref: refs/papers/PDG2024_sum_leptons.pdf (lepton masses)
# Ref: refs/papers/PDG2024_sum_quarks.pdf (quark masses, MS-bar at μ=2 GeV for u,d,s)
# Ref: refs/papers/PDG2024_rev_quark_masses.pdf (c, b at their own scale; t direct)

const SM_LEPTONS = [
    (mass2 = 0.00051100^2, charge2 = 1.0, nc = 1),  # e:  0.51100 MeV
    (mass2 = 0.10566^2,    charge2 = 1.0, nc = 1),  # μ:  105.658 MeV
    (mass2 = 1.77693^2,    charge2 = 1.0, nc = 1),  # τ:  1776.93 MeV
]

const SM_QUARKS = [
    (mass2 = 0.00216^2, charge2 = 4.0/9, nc = 3),  # u:  2.16 MeV (MS-bar, μ=2 GeV)
    (mass2 = 0.00470^2, charge2 = 1.0/9, nc = 3),  # d:  4.70 MeV
    (mass2 = 0.0935^2,  charge2 = 1.0/9, nc = 3),  # s:  93.5 MeV
    (mass2 = 1.2730^2,  charge2 = 4.0/9, nc = 3),  # c:  1.273 GeV (MS-bar, μ=m_c)
    (mass2 = 4.183^2,   charge2 = 1.0/9, nc = 3),  # b:  4.183 GeV (MS-bar, μ=m_b)
    (mass2 = 172.57^2,  charge2 = 4.0/9, nc = 3),  # t:  172.57 GeV (direct meas.)
]

const SM_FERMIONS = vcat(SM_LEPTONS, SM_QUARKS)

"""
    delta_alpha(q2; alpha, fermions) → ComplexF64

Total vacuum polarization Δα(q²) = -Σ_f Π̂_f(q²) summed over active fermions.
Fermions with 4m_f² > |q²| contribute negligibly and are included for completeness.
"""
function delta_alpha(q2::Float64; alpha::Float64 = 1.0 / 137.036,
                     fermions = SM_FERMIONS)::ComplexF64
    result = 0.0 + 0.0im
    for f in fermions
        # Π̂_f = Q_f² × N_c × vacuum_polarization(q², m_f²)
        pi_hat = vacuum_polarization(q2, f.mass2; alpha)
        result -= f.charge2 * f.nc * pi_hat  # Δα = -Π̂
    end
    result
end

"""
    running_alpha(q2; alpha, fermions) → Float64

Running fine-structure constant α(q²) = α / (1 - Δα(q²)).
For spacelike q² < 0 (Euclidean), the result is real.
For timelike q² > 0, returns α / |1 - Δα(q²)| (modulus).
"""
function running_alpha(q2::Float64; alpha::Float64 = 1.0 / 137.036,
                       fermions = SM_FERMIONS)::Float64
    da = delta_alpha(q2; alpha, fermions)
    alpha / abs(1.0 - da)
end

"""
    sigma_improved_ee_mumu(s; alpha) → Float64

Improved Born approximation for σ(e⁺e⁻ → μ⁺μ⁻) including vacuum polarization
resummation and O(α) Schwinger correction. Result in nb.

σ = (4πα²)/(3s) × |1/(1-Π̂(s))|² × (1 + δ_Schwinger)
  = (4π α(s)²)/(3s) × (1 + δ_Schwinger)
"""
function sigma_improved_ee_mumu(s::Float64; alpha::Float64 = 1.0 / 137.036)
    alpha_s = running_alpha(s; alpha)
    sigma_0 = 4.0 * π * alpha_s^2 / (3.0 * s)
    delta = schwinger_correction(; alpha)
    sigma_0 * (1.0 + delta) * 0.3894e6  # GeV⁻² → nb
end
