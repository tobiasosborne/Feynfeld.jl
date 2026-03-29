# c0_analytical.jl — Scalar three-point function C₀
#
# Primary: COLLIER library (Denner, Dittmaier, Hofer, arXiv:1604.06792)
# via ccall to the COLI branch C0_coli function.
# Ref: refs/COLLIER/COLLIER-1.2.8/src/COLI/coli_c0.F, line 35
# Ref: refs/papers/vanOldenborgh1990_ZPhysC46.pdf (FF algorithms)
#
# Special case: C0p0 (all p²=0) — pure logarithms, handled analytically.
# Ref: refs/LoopTools/src/C/C0func.F, lines 401-438
#
# Edge case: Gram-degenerate (Δ₂=0) — COLLIER returns NaN for scalar C₀
# when external momenta are linearly dependent. Falls back to quadgk.
# This is a known limitation of ALL Spence-based algorithms (the 1/√Δ₂
# overall factor diverges). Future: implement the L'Hôpital limit.

const _C0_EPS = 1e-14

# ---- COLLIER library ----
const _COLLIER_LIB = joinpath(@__DIR__, "../../refs/COLLIER/COLLIER-1.2.8/libcollier.so")

let _initialized = Ref(false)
    global function _ensure_collier_init()
        _initialized[] && return
        isfile(_COLLIER_LIB) || return
        ccall((:__collier_init_MOD_init_cll, _COLLIER_LIB), Cvoid,
              (Ref{Int32}, Ref{Int32}, Ptr{UInt8}, Ref{Int32}, Csize_t),
              Ref{Int32}(4), Ref{Int32}(4), "output", Ref{Int32}(0), 6)
        _initialized[] = true
    end
end

# C0_coli(p12,p23,p13, m12,m22,m32) — returns ComplexF64
# Convention mapping: Feynfeld prop 0,1,2 = COLLIER prop 1,2,3 (index shift)
# p10 (props 0-1) → p12_coli (props 1-2), same argument position. No reorder.
function _C0_collier(p10::Float64, p12::Float64, p20::Float64,
                     m02::Float64, m12::Float64, m22::Float64)::ComplexF64
    _ensure_collier_init()
    ccall((:c0_coli_, _COLLIER_LIB), ComplexF64,
          (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
           Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
          ComplexF64(p10), ComplexF64(p12), ComplexF64(p20),
          ComplexF64(m02), ComplexF64(m12), ComplexF64(m22))
end

# ---- C0p0: all momenta² = 0 ----
# Ref: refs/LoopTools/src/C/C0func.F, lines 401-438
function _C0p0(m1::Float64, m2::Float64, m3::Float64)::ComplexF64
    dm12, dm13, dm23 = m1 - m2, m1 - m3, m2 - m3
    eps = _C0_EPS * max(m1, m2, m3, 1.0)
    if abs(dm23) < eps
        if abs(dm13) < eps
            m1 == 0.0 && error("C0(0,0,0, 0,0,0) is IR divergent")
            return complex(-0.5 / m1)
        end
        return complex((dm13 - m1 * log(m1 / m3)) / dm13^2)
    end
    abs(dm12) < eps &&
        return complex((-dm23 + m3 * log(m2 / m3)) / dm23^2)
    abs(dm13) < eps &&
        return complex((dm23 - m2 * log(m2 / m3)) / dm23^2)
    complex(m3 * log(m1 / m3) / (dm13 * dm23) -
            m2 * log(m1 / m2) / (dm12 * dm23))
end

# ---- Dispatcher ----
function _C0_analytical(p10::Float64, p12::Float64, p20::Float64,
                        m02::Float64, m12::Float64, m22::Float64)::ComplexF64
    # All momenta zero → analytical C0p0
    z = (abs(p10) < _C0_EPS) + (abs(p12) < _C0_EPS) + (abs(p20) < _C0_EPS)
    z == 3 && return _C0p0(m02, m12, m22)

    # COLLIER (handles all non-degenerate kinematics)
    if isfile(_COLLIER_LIB)
        result = _C0_collier(p10, p12, p20, m02, m12, m22)
        if isfinite(real(result)) && isfinite(imag(result))
            return result
        end
    end

    # Fallback: quadgk (Gram-degenerate Δ₂≈0, or COLLIER unavailable)
    _C0_quadgk(p10, p12, p20, m02, m12, m22)
end
