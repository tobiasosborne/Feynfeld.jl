# d0_collier.jl — Scalar four-point function D₀
#
# Primary: COLLIER library (Denner, Dittmaier, Hofer, arXiv:1604.06792)
# via ccall to the COLI branch D0_coli function.
# Ref: refs/COLLIER/COLLIER-1.2.8/src/COLI/coli_d0.F, line 37
#
# No quadgk fallback: triple-nested closures cause Julia JIT explosion.
# UV-finite (superficial degree = -4). No μ² dependence.

# D0_coli(p12,p23,p34,p14,p13,p24, m12,m22,m32,m42) → ComplexF64
# Convention mapping: Feynfeld prop 0,1,2,3 = COLLIER vertex 1,2,3,4 (index shift)
# p10 (props 0-1) → p12_coli (verts 1-2), etc. Same argument order, no reorder.
# Ref: refs/COLLIER/COLLIER-1.2.8/src/COLI/coli_d0.F, lines 37-60
function _D0_collier(p10::Float64, p12::Float64, p23::Float64,
                     p30::Float64, p20::Float64, p13::Float64,
                     m02::Float64, m12::Float64, m22::Float64,
                     m32::Float64)::ComplexF64
    _ensure_collier_init()
    ccall((:d0_coli_, _COLLIER_LIB), ComplexF64,
          (Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
           Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64},
           Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}),
          ComplexF64(p10), ComplexF64(p12), ComplexF64(p23),
          ComplexF64(p30), ComplexF64(p20), ComplexF64(p13),
          ComplexF64(m02), ComplexF64(m12), ComplexF64(m22), ComplexF64(m32))
end

# Dispatcher: COLLIER only (no quadgk fallback — triple-nested closures
# cause a JIT compilation explosion in Julia, and take O(minutes) at runtime).
function _D0_evaluate(p10::Float64, p12::Float64, p23::Float64,
                      p30::Float64, p20::Float64, p13::Float64,
                      m02::Float64, m12::Float64, m22::Float64,
                      m32::Float64)::ComplexF64
    isfile(_COLLIER_LIB) ||
        error("D₀ requires COLLIER library at $_COLLIER_LIB — see HANDOFF.md for build instructions")
    result = _D0_collier(p10, p12, p23, p30, p20, p13,
                         m02, m12, m22, m32)
    (isfinite(real(result)) && isfinite(imag(result))) ||
        error("COLLIER D₀ returned non-finite: $result")
    result
end
