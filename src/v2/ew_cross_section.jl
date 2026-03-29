# Tree-level e+e- → W+W- total cross-section (massless electron limit).
#
# 3 diagrams: s-channel γ, s-channel Z, t-channel ν_e.
# (Higgs diagram vanishes for m_e = 0.)
#
# Ref: refs/FeynCalc/FeynCalc/Examples/EW/Tree/Mathematica/AnelEl-WW.m
#      lines 193-198: verified "CORRECT" against Grozin formula
# Secondary: Grozin, "Using REDUCE in HEP", Ch. 5.4 (not locally stored)
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.1)-(11.2)
#      (process definition and kinematics for e+e-→W+W-)
#
# The formula uses s in units of M_W² (dimensionless):
#   σ = Part_log(s) + Part_sqrt(s)    [in units of M_W⁻²]
#
# Part_log has a logarithm:
#   L(s) = ln[(s-2-√((s-4)s)) / (s-2+√((s-4)s))]
#
# Part_sqrt has √(s-4) (W velocity factor).
#
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Eqs. (11.1)-(11.2)
# "e⁺(p₁) + e⁻(p₂) → W⁺(k₁) + W⁻(k₂)"
# "s = (p₁+p₂)², t = (p₁-k₁)², s+t+u = 2m_e²+2M_W²"

# ──────────────────────────────────────────────────────────────────
# REFERENCE IMPLEMENTATION — validates physics formula directly.
# The pipeline should reproduce this result via
# Model → Rules → Diagrams → Algebra → Integrals → Evaluate.
# Do NOT delete: used for cross-validation against pipeline results.
# ──────────────────────────────────────────────────────────────────

"""
    sigma_ee_ww(s_gev2; M_W, sin2_W, alpha) → Float64

Tree-level total cross-section σ(e⁺e⁻ → W⁺W⁻) in picobarns.
Massless electron limit (m_e = 0), no Higgs diagram.

`s_gev2`: center-of-mass energy squared in GeV².
"""
function sigma_ee_ww(s_gev2::Float64;
                     M_W::Float64 = EW_M_W,
                     sin2_W::Float64 = EW_SIN2_W,
                     alpha::Float64 = EW_ALPHA)::Float64
    s = s_gev2 / M_W^2          # dimensionless, in M_W² units
    s <= 4.0 && return 0.0       # below threshold 2M_W
    _sigma_grozin(s, sin2_W, alpha) / M_W^2 * 3.894e8  # M_W⁻² → GeV⁻² → pb
end

# Grozin formula: total cross-section in M_W⁻² units.
# Ref: FeynCalc AnelEl-WW.m, lines 193-198
# Ref: Grozin, "Using REDUCE in High Energy Physics", Chapter 5.4
function _sigma_grozin(s::Float64, sin2_W::Float64, alpha::Float64)::Float64
    cos2_W = 1.0 - sin2_W
    sin4_W = sin2_W^2
    a2 = alpha^2

    # W velocity
    sq = sqrt((s - 4.0) * s)

    # Logarithm: L = ln[(s-2-√((s-4)s)) / (s-2+√((s-4)s))]
    L = log((s - 2.0 - sq) / (s - 2.0 + sq))

    # Z-pole factor: (1 - s cos²θ_W) = (M_W² - s cos²θ_W M_W²)/M_W²
    # = (M_W² - s_phys cos²θ_W)/M_W² ∝ (s_phys - M_Z²) at the Z pole
    zpole = 1.0 - s * cos2_W   # negative for s > 1/cos²θ_W (above Z)

    # Part 1: logarithmic term
    poly_log = 24.0 * s * (4.0 + s + s^2) -
               24.0 * (4.0 + 10.0*s + 2.0*s^2 + s^3) * sin2_W
    part_log = π * L * a2 * poly_log / (96.0 * s^3 * zpole * sin4_W)

    # Part 2: polynomial × √(s-4) term
    poly_sq = -3.0 * s * (32.0 - 20.0*s + 21.0*s^2) +
               16.0 * (3.0 + 8.0*s) * (2.0 + s^2) * sin2_W -
               4.0 * (96.0 + 160.0*s + 8.0*s^2 + 15.0*s^3) * sin4_W
    part_sq = π * sqrt(s - 4.0) * a2 * poly_sq /
              (96.0 * s^2.5 * zpole^2 * sin4_W)

    part_log + part_sq
end

"""
    dsigma_dt_ee_ww(s_gev2, t_gev2; M_W, sin2_W, alpha) → Float64

Tree-level differential cross-section dσ/dt for e⁺e⁻ → W⁺W⁻ in pb/GeV².
(Not yet implemented — future extension.)
"""
function dsigma_dt_ee_ww(s_gev2::Float64, t_gev2::Float64; kwargs...)
    error("dsigma_dt_ee_ww not yet implemented; use sigma_ee_ww for total")
end
