# QED vertex correction: anomalous magnetic moment F₂(0).
#
# Ground truth: F₂(0) = α/(2π)  (Schwinger, Phys. Rev. 73 (1948) 416)
# Cross-check: refs/FeynCalc/FeynCalc/Examples/QED/OneLoop/Mathematica/El-GaEl.m
#              (line 210: knownResult = AlphaFS/(2*Pi))
#
# Ref: refs/papers/(Frontiers_in_Physics)...Peskin...Schroeder...(1995).djvu
#      Chapter 6, Eqs. 6.47-6.56 (.djvu unreadable by tools; verified via El-GaEl.m)
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Appendix C, Eqs. (C.1)-(C.39)
#      (general vertex form factors; our F₂ formula is the QED on-shell specialization)
#
# After Dirac algebra, on-shell projection, and D-dimensional loop integration:
#   F₂(q²) = (α/π) m² ∫₀¹ dz ∫₀^{1-z} dx  z(1-z) / Δ(x,z)
#   Δ(x,z) = m²(1-z)² - q² x(1-x-z) + z λ²
#
# At q²=0 the integrand depends only on z (not x), so the x-integral gives (1-z):
#   F₂(0) = (α/π) m² ∫₀¹ dz  z(1-z)² / [m²(1-z)² + z λ²]
#
# For λ→0 (IR-finite): F₂(0) = (α/π) × ½ = α/(2π)

using QuadGK: quadgk

"""
    vertex_f2_zero(m2, lambda2; alpha) → Float64

Anomalous magnetic moment F₂(q²=0) from QED vertex correction at one loop.

`m2`: electron mass squared, `lambda2`: photon mass squared (IR regulator),
`alpha`: fine-structure constant (default 1/137.036).
"""
function vertex_f2_zero(m2::Float64, lambda2::Float64;
                        alpha::Float64 = 1.0 / 137.036)::Float64
    val, _ = quadgk(0.0, 1.0; rtol = 1e-12) do z
        m2 * z * (1.0 - z)^2 / (m2 * (1.0 - z)^2 + z * lambda2)
    end
    (alpha / π) * val
end

"""
    vertex_f2(q2, m2, lambda2; alpha) → Float64

Dirac form factor F₂(q²) from QED vertex correction at one loop.

General q² version (2D Feynman parameter integral). For q²=0, prefer
`vertex_f2_zero` which reduces to a 1D integral.

Valid for q² ≤ 0 (spacelike) and 0 < q² < 4m² (below threshold).
Above threshold the iε prescription is not implemented.
"""
function vertex_f2(q2::Float64, m2::Float64, lambda2::Float64;
                   alpha::Float64 = 1.0 / 137.036)::Float64
    val, _ = quadgk(0.0, 1.0; rtol = 1e-10) do z
        inner, _ = quadgk(0.0, 1.0 - z; rtol = 1e-10) do x
            delta = m2 * (1.0 - z)^2 - q2 * x * (1.0 - x - z) + z * lambda2
            m2 * z * (1.0 - z) / delta
        end
        inner
    end
    (alpha / π) * val
end
