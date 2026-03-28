# Polarization sum for gauge bosons.
#
# Ref: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID1
# Ref: Peskin & Schroeder, discussion around Eq. (5.75)
#
# For a massless gauge boson (photon) in Feynman gauge:
#   Σ_λ ε^μ(k,λ) ε^ν*(k,λ) = -g^{μν}
#
# This is the only form needed for tree-level QED. Massive boson and
# axial/light-cone gauge forms are deferred to later spirals.

"""
    polarization_sum(mu, nu)

Feynman gauge polarization sum: Σ_λ ε^μ ε^{ν*} = -g^{μν}.
"""
polarization_sum(mu::LorentzIndex, nu::LorentzIndex) = -alg(pair(mu, nu))
