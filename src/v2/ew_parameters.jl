# Electroweak Standard Model parameters and derived couplings.
#
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Eqs. (10.11), (10.22a), (10.63)
#      "α⁻¹ = 137.036", "sin²θ_W = 0.22348", "M_W = 80.360 GeV"
# Ref: refs/papers/PDG2024_list_z_boson.pdf
#      "M_Z = 91.1880 ± 0.0020 GeV"

# ---- Physical constants (PDG 2024 on-shell scheme) ----

const EW_M_W  = 80.360     # W boson mass [GeV]
const EW_M_Z  = 91.1880    # Z boson mass [GeV]
const EW_SIN2_W = 0.22348  # sin²θ_W (on-shell)
const EW_COS2_W = 1.0 - EW_SIN2_W
const EW_SIN_W = sqrt(EW_SIN2_W)
const EW_COS_W = sqrt(EW_COS2_W)
const EW_ALPHA = 1.0 / 137.036

# ---- Derived couplings ----

# Z-fermion vector/axial couplings for electron (T₃ = -1/2, Q = -1):
# Ref: refs/papers/PDG2024_rev_standard_model.pdf, Table 10.3
# "g_V^e = -1/2 + 2sin²θ_W", "g_A^e = -1/2"
const EW_GV_E = -0.5 + 2.0 * EW_SIN2_W
const EW_GA_E = -0.5

# Z-fermion couplings in left/right form:
# g_L = g_V + g_A = -1 + 2sin²θ_W,  g_R = g_V - g_A = 2sin²θ_W
const EW_GL_E = EW_GV_E + EW_GA_E
const EW_GR_E = EW_GV_E - EW_GA_E
