# Feynfeld.jl — Passarino-Veltman scalar integral symbols
#
# PaVe{N} represents an N-point scalar or tensor coefficient function:
#   PaVe{1} = A-functions (1-point: A0, A00)
#   PaVe{2} = B-functions (2-point: B0, B1, B00, B11)
#   PaVe{3} = C-functions (3-point: C0, C1, C2, ...)
#   PaVe{4} = D-functions (4-point: D0, D1, ...)
#
# Convention: PaVe{N}(indices, invariants, masses) where
#   indices: sorted tensor coefficient labels (Int[])
#   invariants: N(N-1)/2 kinematic invariants (Denner convention)
#   masses: N squared masses
#
# Ref: FeynCalc LoopIntegrals/PaVe.m, Denner (1993)

export PaVe, A0, A00, B0, B1, B00, B11, C0, D0

"""
    PaVe{N}(indices, invariants, masses)

Passarino-Veltman coefficient function for an N-point integral.
Indices are automatically sorted (symmetry of tensor decomposition).

# Examples
```julia
A0(:m²)                              # PaVe{1}
B0(:pp, :m1², :m2²)                  # PaVe{2}
PaVe{3}([1,2], [:p10,:p12,:p20], [:m1²,:m2²,:m3²])  # C12
```
"""
struct PaVe{N} <: FeynExpr
    indices::Vector{Int}
    invariants::Vector{Any}
    masses::Vector{Any}
    function PaVe{N}(indices::Vector{Int}, invariants::Vector{Any},
                     masses::Vector{Any}) where {N}
        N >= 1 || error("PaVe{N}: N must be >= 1, got $N")
        # 1-point functions have no tensor indices with value >= 1
        if N == 1 && any(i -> i >= 1, indices)
            error("PaVe{1}: tensorial 1-point functions with index >= 1 do not exist")
        end
        n_inv = N * (N - 1) ÷ 2
        length(invariants) == n_inv || error(
            "PaVe{$N}: expected $n_inv invariants, got $(length(invariants))")
        length(masses) == N || error(
            "PaVe{$N}: expected $N masses, got $(length(masses))")
        new{N}(sort(indices), invariants, masses)
    end
end

function Base.:(==)(a::PaVe{N}, b::PaVe{N}) where {N}
    a.indices == b.indices && a.invariants == b.invariants && a.masses == b.masses
end
Base.:(==)(::PaVe{N}, ::PaVe{M}) where {N,M} = false

Base.hash(p::PaVe{N}, h::UInt) where {N} =
    hash(p.masses, hash(p.invariants, hash(p.indices, hash(N, hash(:PaVe, h)))))

# ── Named constructors ──────────────────────────────────────────────

# 1-point (A-functions)
"""A0(m²) — scalar 1-point function."""
A0(m) = PaVe{1}(Int[], Any[], Any[m])
"""A00(m²) — tensor 1-point function (g^{μν} coefficient)."""
A00(m) = PaVe{1}([0, 0], Any[], Any[m])

# 2-point (B-functions)
"""B0(p², m0², m1²) — scalar 2-point function."""
B0(p, m0, m1) = PaVe{2}(Int[], Any[p], Any[m0, m1])
"""B1(p², m0², m1²) — tensor 2-point, p^μ coefficient."""
B1(p, m0, m1) = PaVe{2}([1], Any[p], Any[m0, m1])
"""B00(p², m0², m1²) — tensor 2-point, g^{μν} coefficient."""
B00(p, m0, m1) = PaVe{2}([0, 0], Any[p], Any[m0, m1])
"""B11(p², m0², m1²) — tensor 2-point, p^μ p^ν coefficient."""
B11(p, m0, m1) = PaVe{2}([1, 1], Any[p], Any[m0, m1])

# 3-point (C-functions)
"""C0(p10², p12², p20², m0², m1², m2²) — scalar 3-point function."""
C0(p10, p12, p20, m0, m1, m2) =
    PaVe{3}(Int[], Any[p10, p12, p20], Any[m0, m1, m2])

# 4-point (D-functions)
"""D0(p10², p12², p23², p30², p20², p13², m0², m1², m2², m3²) — scalar 4-point."""
function D0(p10, p12, p23, p30, p20, p13, m0, m1, m2, m3)
    PaVe{4}(Int[], Any[p10, p12, p23, p30, p20, p13], Any[m0, m1, m2, m3])
end
