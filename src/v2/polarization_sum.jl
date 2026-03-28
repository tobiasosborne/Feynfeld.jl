# Polarization sum for gauge bosons.
#
# Ref: refs/FeynCalc/Tests/Feynman/PolarizationSum.test, ID1
# Ref: Peskin & Schroeder, discussion around Eq. (5.75)
#
# Feynman gauge: Σ_λ ε^μ ε^{ν*} = -g^{μν}
# Axial gauge:   Σ_λ ε^μ ε^{ν*} = -g^{μν} + (k^μ n^ν + n^μ k^ν)/(k·n)
#   for massless k and massless reference momentum n (n²=0, k·n≠0).
#
# The axial gauge form sums over physical (transverse) polarizations only.
# Required for QCD where individual diagram contributions are not gauge-invariant.

"""
    polarization_sum(mu, nu)

Feynman gauge polarization sum: Σ_λ ε^μ ε^{ν*} = -g^{μν}.
"""
polarization_sum(mu::LorentzIndex, nu::LorentzIndex) = -alg(pair(mu, nu))

"""
    polarization_sum(mu, nu, k, n)

Axial gauge polarization sum for massless gauge boson with momentum `k`
and massless reference momentum `n` (n²=0):
  Σ_λ ε^μ ε^{ν*} = -g^{μν} + (k^μ n^ν + n^μ k^ν)/(k·n)

Ref: FeynCalc DoPolarizationSums with reference vector.
"""
function polarization_sum(mu::LorentzIndex, nu::LorentzIndex,
                          k::Momentum, n::Momentum; ctx::SPContext=CURRENT_SP[])
    kn = _lookup_sp(k, n, ctx)
    kn === nothing && error("polarization_sum: need k·n in SPContext for $(k.name)·$(n.name)")
    -alg(pair(mu, nu)) +
        (1 // kn) * (alg(pair(mu, k)) * alg(pair(nu, n)) +
                     alg(pair(mu, n)) * alg(pair(nu, k)))
end

"""
    polarization_sum_massive(mu, nu, k, M2; ctx)

Massive vector boson polarization sum:
  Σ_λ ε^μ ε^{ν*} = -g^{μν} + k^μ k^ν / M²

Ref: Peskin & Schroeder, Eq. (5.75) generalized to massive case.
Used for external W± and Z bosons.
"""
function polarization_sum_massive(mu::LorentzIndex, nu::LorentzIndex,
                                  k::Momentum, M2::Number)
    -alg(pair(mu, nu)) + (1 // M2) * alg(pair(mu, k)) * alg(pair(nu, k))
end

function _lookup_sp(a::Momentum, b::Momentum, ctx::SPContext)
    key = a.name <= b.name ? (a.name, b.name) : (b.name, a.name)
    get(ctx.values, key, nothing)
end
