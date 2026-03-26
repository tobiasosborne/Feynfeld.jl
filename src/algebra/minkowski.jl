# Feynfeld.jl — Minkowski spacetime registry (TensorGR bridge)
#
# Sets up a TensorGR registry with a Minkowski manifold, flat metric η,
# and Lorentz signature. This is the bridge between Feynfeld's QFT algebra
# and TensorGR's index contraction/canonicalization engine.

export minkowski_registry

"""
    minkowski_registry(; dim=4, signature=:mostly_plus)

Create a TensorGR registry pre-configured for Minkowski spacetime.

The metric η has Lorentz signature (default: mostly-plus, i.e. −+++ convention).
For dimensional regularisation, pass `dim=:D`.

Returns a `TensorGR.TensorRegistry`.
"""
function minkowski_registry(; dim::Union{Int,Symbol}=4, signature=:mostly_plus)
    reg = TensorGR.TensorRegistry()

    # Index alphabet for Lorentz indices
    lorentz_indices = [:μ, :ν, :ρ, :σ, :τ, :κ, :λ, :α, :β, :γ, :δ, :ε]

    # Register Minkowski manifold
    mp = TensorGR.ManifoldProperties(:M4, dim, :η, nothing, lorentz_indices)
    TensorGR.register_manifold!(reg, mp)

    # Register the Minkowski metric η
    sig = if dim isa Int
        if signature === :mostly_plus
            TensorGR.lorentzian(dim)
        elseif signature === :mostly_minus
            TensorGR.MetricSignature(vcat([1], fill(-1, dim - 1)))
        else
            error("Unknown signature: $signature. Use :mostly_plus or :mostly_minus")
        end
    else
        nothing  # symbolic dim: no concrete signature
    end
    TensorGR.define_metric!(reg, :η; manifold=:M4, signature=sig)

    # Minkowski is flat: Riemann = Ricci = RicciScalar = Christoffel = 0
    TensorGR.set_flat!(reg, :η)

    reg
end
