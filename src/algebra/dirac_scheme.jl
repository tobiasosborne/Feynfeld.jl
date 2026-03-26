# Feynfeld.jl — Dirac gamma5 scheme system
#
# Supported schemes:
#   :NDR   — Naive Dimensional Regularisation (default)
#            gamma5 anticommutes with all gamma^mu in D dimensions
#   :BMHV  — Breitenlohner-Maison / 't Hooft-Veltman
#            gamma5 anticommutes with 4D gammas only, not evanescent
#   :Larin — Larin scheme
#            gamma5 replaced by (i/6) ε^{μνρσ} γ_μ γ_ν γ_ρ in traces
#
# For Phase 1b, only NDR is fully implemented. BMHV and Larin are
# registered but defer to NDR for most operations.
#
# Ref: FeynCalc Shared/SharedTools.m $BreitMaison, $Larin

export DiracScheme, NDR, BMHV_SCHEME, LARIN, current_scheme, set_scheme!, with_scheme

"""Dirac gamma5 regularisation scheme."""
@enum DiracScheme NDR=1 BMHV_SCHEME=2 LARIN=3

"""Current global scheme (thread-local in future)."""
const _DIRAC_SCHEME = Ref{DiracScheme}(NDR)

"""Get the current Dirac scheme."""
current_scheme() = _DIRAC_SCHEME[]

"""Set the Dirac scheme."""
function set_scheme!(s::DiracScheme)
    _DIRAC_SCHEME[] = s
    s
end

"""
    with_scheme(f, scheme)

Execute `f()` with a temporary Dirac scheme, restoring the original after.
"""
function with_scheme(f, scheme::DiracScheme)
    old = _DIRAC_SCHEME[]
    _DIRAC_SCHEME[] = scheme
    try
        f()
    finally
        _DIRAC_SCHEME[] = old
    end
end
