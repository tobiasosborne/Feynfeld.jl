# c0_analytical.jl — Analytical scalar three-point function C₀
#
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Section 5, Appendix B
#   "'t Hooft, Veltman: Scalar One-Loop Integrals (1979)"
#   C₀ definition: Eq. (5.1), Feynman param: Eq. (5.2)
#   Spence decomposition: Eqs. (B.1)-(B.3), η function: Eq. (2.4)
# Ref: refs/LoopTools/src/C/C0func.F (Denner adaptation)
#   Dispatch: lines 108-119, C0p0: 401-438, C0p1: 334-397,
#   C0p2: 266-329, C0p3: 168-261
#
# Dispatch by number of zero external momenta²:
#   C0p0 (z=3): pure logarithms, 5 mass subcases
#   C0p1 (z=2): Li₂ with roots from one quadratic
#   C0p2 (z=1): Li₂ with different variable setup
#   C0p3 (z=0): full Li₂ + η-function corrections
# For z ∈ {0,1,2}: C₀ = Σ_{cyclic perms} channel(perm)

using PolyLog: li2

# LoopTools spence(0,x,0) = Li₂(x). Internally: Li2series(1-x) = Li₂(1-(1-x)) = Li₂(x).
# Ref: refs/LoopTools/src/util/auxCD.F, lines 262-361
_sp(x::Number) = li2(x)

const _C0_EPS  = 1e-14  # threshold for "zero momentum²"
const _C0_IEPS = 1e-50  # infinitesimal imaginary part (iε)

# ---- C0p0: all momenta² = 0 ----
# Ref: refs/LoopTools/src/C/C0func.F, lines 401-438
# "C0p0 = m3·ln(m1/m3)/(m13·m23) - m2·ln(m1/m2)/(m12·m23)  [all distinct]"
# All arguments are squared masses.
function _C0p0(m1::Float64, m2::Float64, m3::Float64)::ComplexF64
    dm12, dm13, dm23 = m1 - m2, m1 - m3, m2 - m3
    eps = _C0_EPS * max(m1, m2, m3, 1.0)
    if abs(dm23) < eps
        if abs(dm13) < eps  # all equal: C₀ = -1/(2m²)
            m1 == 0.0 && error("C0(0,0,0, 0,0,0) is IR divergent")
            return complex(-0.5 / m1)
        end
        return complex((dm13 - m1 * log(m1 / m3)) / dm13^2)  # m2 ≈ m3
    end
    abs(dm12) < eps &&
        return complex((-dm23 + m3 * log(m2 / m3)) / dm23^2)  # m1 ≈ m2
    abs(dm13) < eps &&
        return complex((dm23 - m2 * log(m2 / m3)) / dm23^2)   # m1 ≈ m3
    # All distinct
    complex(m3 * log(m1 / m3) / (dm13 * dm23) -
            m2 * log(m1 / m2) / (dm12 * dm23))
end

# ---- C0p1: two momenta² = 0 ----
# Ref: refs/LoopTools/src/C/C0func.F, lines 334-397
# px1 is the single nonzero momentum². Returns 0 if px1 ≈ 0.
function _C0p1_channel(px1::Float64, px2::Float64, px3::Float64,
                       mx1::Float64, mx2::Float64, mx3::Float64)::ComplexF64
    abs(px1) < _C0_EPS && return 0.0 + 0.0im
    dm12, dm13, dm23 = mx1 - mx2, mx1 - mx3, mx2 - mx3

    # Piece A: logarithmic contribution (only if dm13 ≠ 0)
    piece_a = 0.0 + 0.0im
    if abs(dm13) > _C0_EPS * max(abs(mx1), abs(mx3), 1.0)
        ya1 = complex(dm23 - px1)
        ya2 = complex(dm23)
        c = dm23 + px1 * mx3 / dm13
        ya3 = complex(c, -copysign(abs(c), px1 / dm13) * _C0_IEPS)
        piece_a = _sp(ya1 / ya3) - _sp(ya2 / ya3)
    end

    # Piece B: Spence from quadratic roots
    y1 = complex(-2.0 * px1 * dm23)
    y2 = complex(-2.0 * px1 * (dm23 - px1))
    c = px1 * (px1 - dm13 - dm23)
    b = px1 * sqrt(complex((px1 - dm12)^2 - 4.0 * px1 * mx2))
    y3, y4 = c - b, c + b
    # Stability: y3·y4 = 4px1²(px1·mx3 + dm13·dm23)
    n34 = complex(4.0 * px1^2 * (px1 * mx3 + dm13 * dm23))
    if abs(y3) < abs(y4); y3 = n34 / y4
    else                 ; y4 = n34 / y3; end
    y3 -= abs(y3) * complex(0.0, _C0_IEPS)
    y4 += abs(y4) * complex(0.0, _C0_IEPS)

    (piece_a + _sp(y1/y3) + _sp(y1/y4) - _sp(y2/y3) - _sp(y2/y4)) / px1
end

# ---- C0p2: one momentum² = 0 ----
# Ref: refs/LoopTools/src/C/C0func.F, lines 266-329
# px1 nonzero; one of px2, px3 is zero.
function _C0p2_channel(px1::Float64, px2::Float64, px3::Float64,
                       mx1::Float64, mx2::Float64, mx3::Float64)::ComplexF64
    abs(px1) < _C0_EPS && return 0.0 + 0.0im
    dm12, dm13, dm23 = mx1 - mx2, mx1 - mx3, mx2 - mx3

    if abs(px3) < _C0_EPS
        a = px1 - px2
        y1 = complex(-2.0 * px1 * (dm13 - px1 + px2))
        y2 = complex(-2.0 * px1 * dm13)
    else
        a = px3 - px1
        y1 = complex(-2.0 * px1 * dm23)
        y2 = complex(-2.0 * px1 * (dm23 + px3 - px1))
    end
    abs(a) < _C0_EPS && return 0.0 + 0.0im

    c = px1 * (px1 - px2 - px3 - dm13 - dm23) - dm12 * (px2 - px3)
    b = a * sqrt(complex((px1 - dm12)^2 - 4.0 * px1 * mx2))
    y3, y4 = c + b, c - b
    # Stability: known product (C0func.F lines 312-320)
    n34 = complex(4.0 * px1 * (
        px1 * ((px1 - px2 - px3) * mx3 + px2 * px3 + dm13 * dm23) +
        px2 * ((px2 - px3 - px1) * mx1 + dm12 * dm13) +
        px3 * ((px3 - px1 - px2) * mx2 - dm12 * dm23)))
    if abs(y3) < abs(y4); y3 = n34 / y4
    else                 ; y4 = n34 / y3; end
    sgn = copysign(1.0, a / px1)
    y3 += complex(0.0, sgn * abs(y3) * _C0_IEPS)
    y4 -= complex(0.0, sgn * abs(y4) * _C0_IEPS)

    (_sp(y2/y3) + _sp(y2/y4) - _sp(y1/y3) - _sp(y1/y4)) / a
end

# ---- C0p3: all momenta² nonzero ----
# Ref: refs/LoopTools/src/C/C0func.F, lines 168-261
# Spence decomposition (Denner adaptation).
# NOTE: The Denner η-function corrections (C0func.F lines 234-256) are omitted.
# C0func.F is a backup double-checker for the FF library (ffxc0), and the η
# corrections are numerically fragile when Kallen(p1,p2,p3) < 0 (complex a).
# The η-free formula is exact for spacelike/below-threshold kinematics.
# For above-threshold timelike momenta, the dispatcher falls back to quadgk.
function _C0p3_channel(px1::Float64, px2::Float64, px3::Float64,
                       mx1::Float64, mx2::Float64, mx3::Float64)::ComplexF64
    dm12, dm13, dm23 = mx1 - mx2, mx1 - mx3, mx2 - mx3
    a2 = (px1 - px2 - px3)^2 - 4.0 * px2 * px3  # Kallen(px1,px2,px3)
    a = sqrt(complex(a2))
    n = 0.5 / px1
    cv = (px1 * (px1 - px2 - px3 - dm13 - dm23) - dm12 * (px2 - px3)) / a

    # Products for stability (C0func.F lines 202-209)
    n123 = px1 * (px2 * px3 + dm13 * dm23) + dm12 * (dm13 * px2 - dm23 * px3)
    pp1 = px1 * (px1 - px2 - px3)
    pp2 = px2 * (px1 - px2 + px3)
    pp3 = px3 * (px1 + px2 - px3)
    n1 = n123 - dm23 * pp1 - dm12 * pp2
    n2 = n123 - dm13 * pp1 + dm12 * pp3
    n3 = n123 + mx3 * pp1 - mx1 * pp2 - mx2 * pp3

    # Four roots (C0func.F lines 211-229)
    y1  = n * (cv + (px1 - dm12))
    y1a = n * (cv - (px1 - dm12))
    if abs(y1) < abs(y1a); y1 = n1 / (a^2 * px1 * y1a); end
    y2  = n * (cv - (px1 + dm12))
    y2a = n * (cv + (px1 + dm12))
    if abs(y2) < abs(y2a); y2 = n2 / (a^2 * px1 * y2a); end
    bv = sqrt(complex((px1 - dm12)^2 - 4.0 * px1 * mx2))
    y3 = n * (cv + bv)
    y4 = n * (cv - bv)
    if abs(y3) < abs(y4); y3 = n3 / (a^2 * px1 * y4)
    else                 ; y4 = n3 / (a^2 * px1 * y3); end

    s = real(a * bv)  # iε sign
    y3 += complex(0.0, copysign(abs(y3), s) * _C0_IEPS)
    y4 -= complex(0.0, copysign(abs(y4), s) * _C0_IEPS)

    (_sp(y2/y3) + _sp(y2/y4) - _sp(y1/y3) - _sp(y1/y4)) / a
end

# ---- Threshold check ----
# The η-free Denner formula is exact below threshold. Above threshold,
# the Spence arguments can cross branch cuts and the formula needs η
# corrections which are numerically fragile. Fall back to quadgk instead.
function _C0_above_threshold(p10::Float64, p12::Float64, p20::Float64,
                             m02::Float64, m12::Float64, m22::Float64)::Bool
    # Threshold for each channel: p > (√m_i + √m_j)²
    t01 = (sqrt(max(m02, 0.0)) + sqrt(max(m12, 0.0)))^2
    t02 = (sqrt(max(m02, 0.0)) + sqrt(max(m22, 0.0)))^2
    t12 = (sqrt(max(m12, 0.0)) + sqrt(max(m22, 0.0)))^2
    # Use 0.99× to catch "at threshold" cases too
    p10 > 0.99 * t01 || p20 > 0.99 * t02 || p12 > 0.99 * t12
end

# ---- Dispatcher ----
# Convention mapping (Feynfeld ↔ LoopTools):
#   p10 ↔ P(1), p20 ↔ P(2), p12 ↔ P(3)
#   m02 ↔ M(1), m12 ↔ M(2), m22 ↔ M(3)
# Cyclic permutations: (1,2,3) → (2,3,1) → (3,1,2)
function _C0_analytical(p10::Float64, p12::Float64, p20::Float64,
                        m02::Float64, m12::Float64, m22::Float64)::ComplexF64
    z = (abs(p10) < _C0_EPS) + (abs(p12) < _C0_EPS) + (abs(p20) < _C0_EPS)
    z == 3 && return _C0p0(m02, m12, m22)
    # Above threshold: η-free formula unreliable → quadgk
    if _C0_above_threshold(p10, p12, p20, m02, m12, m22)
        return _C0_quadgk(p10, p12, p20, m02, m12, m22)
    end
    f = z == 2 ? _C0p1_channel : z == 1 ? _C0p2_channel : _C0p3_channel
    result = f(p10, p20, p12, m02, m12, m22) +  # perm (1,2,3)
             f(p20, p12, p10, m12, m22, m02) +  # perm (2,3,1)
             f(p12, p10, p20, m22, m02, m12)    # perm (3,1,2)
    # Fall back for degenerate cases (a≈0 in C0p2/C0p3)
    if !isfinite(result) || (result == 0.0 + 0.0im && z < 3)
        return _C0_quadgk(p10, p12, p20, m02, m12, m22)
    end
    result
end
