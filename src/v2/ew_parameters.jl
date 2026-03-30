# Electroweak Standard Model parameters and derived couplings.
#
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eqs. (10.11), (10.22a), (10.63)
#      "α⁻¹ = 137.035999178(8)"                     [Eq. (10.11)]
#      "sin²θ_W = 1 - M_W²/M_Z²"                    [Eq. (10.22a), on-shell definition]
#      "sin²θ_W = 0.22348 ± 0.00006"                 [Table 10.2, on-shell]
# Ref: refs/papers/PDG2024_list_z_boson.pdf
#      "M_Z = 91.1880 ± 0.0020 GeV"
# M_W derived: M_W = M_Z √(1 - sin²θ_W) = 91.1880 × √(0.77652) = 80.360 GeV

# ---- Physical constants (PDG 2024 on-shell scheme) ----
# Float64 for numerics, Rational for symbolic algebra (avoids Int64 overflow).

const EW_M_W  = 80.360     # W boson mass [GeV]
const EW_M_Z  = 91.1880    # Z boson mass [GeV]
const EW_SIN2_W = 0.22348  # sin²θ_W (on-shell) — Float64 for numerics
const EW_COS2_W = 1.0 - EW_SIN2_W
const EW_SIN_W = sqrt(EW_SIN2_W)
const EW_COS_W = sqrt(EW_COS2_W)
const EW_ALPHA = 1.0 / 137.036

# Exact Rationals for symbolic algebra (small denominators → no overflow)
const EW_SIN2_W_R = 5587//25000     # = 0.22348 exact
const EW_COS2_W_R = 1//1 - EW_SIN2_W_R

# ---- Derived couplings ----

# Z-fermion vector/axial couplings for electron (T₃ = -1/2, Q = -1):
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Table 10.3
# "g_V^e = -1/2 + 2sin²θ_W", "g_A^e = -1/2"
const EW_GV_E = -0.5 + 2.0 * EW_SIN2_W
const EW_GA_E = -0.5
const EW_GV_E_R = -1//2 + 2 * EW_SIN2_W_R   # Rational version: -663//12500
const EW_GA_E_R = -1//2                        # Rational version

# Z-fermion couplings in left/right form:
# g_L = g_V + g_A = -1 + 2sin²θ_W,  g_R = g_V - g_A = 2sin²θ_W
const EW_GL_E = EW_GV_E + EW_GA_E
const EW_GR_E = EW_GV_E - EW_GA_E
