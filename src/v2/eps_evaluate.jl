# Numerical evaluation of fully-contracted Levi-Civita tensors.
#
# After Lorentz contraction, Eps factors may have all-Momentum slots:
# ε(p1, k2, p2, k1) = ε_{μνρσ} p1^μ k2^ν p2^ρ k1^σ.
# These are parity-odd scalars that need numerical evaluation.
#
# Ref: refs/papers/MertigBohmDenner1991_FeynCalc_CPC64.pdf, Eq. (2.21)
#   "ε^{μνρσ} ε_{μ'ν'ρ'σ'} = -det|g^μ_{μ'} ...|"
#   ⟹ ε(a,b,c,d)² = -det[G] where G_{ij} = a_i · a_j (Gram matrix)
#
# Sign determination: |ε| is determined by the Gram determinant, but the
# sign of ε CANNOT be derived from scalar products alone — it requires an
# orientation convention.
# Ref: refs/papers/Denner1993_FortschrPhys41.pdf, Sec. 11 (implicit)
#
# Convention: ε(p1,p2,k1,k2) > 0 when θ ∈ (0,π) in the CM frame.
# Permutations of the reference order (:p1,:p2,:k1,:k2) acquire parity signs.
# This matches ε = 2|p||p'|√s sinθ > 0 for standard scattering geometry.

using LinearAlgebra: det

# ---- Dispatch-based factor evaluation for evaluate_numeric ----

_eval_factor(f::Pair{Momentum, Momentum}, sp::Dict{Tuple{Symbol,Symbol}, Float64}) =
    sp[_sp_key(f.a.name, f.b.name)]

function _eval_factor(f::Eps, sp::Dict{Tuple{Symbol,Symbol}, Float64})
    _eps_all_momentum(f) || error("evaluate_numeric: Eps with Lorentz indices: $f")
    _evaluate_eps(f, sp)
end

_eval_factor(f, ::Dict{Tuple{Symbol,Symbol}, Float64}) =
    error("evaluate_numeric: unevaluated factor $f")

# ---- Eps evaluation ----

# True if all four Eps slots are Momentum (not LorentzIndex or MomentumSum).
# MomentumSum slots should have been expanded upstream by expand_scalar_product.
_eps_all_momentum(e::Eps) = e.a isa Momentum && e.b isa Momentum &&
                             e.c isa Momentum && e.d isa Momentum

# Evaluate ε(a,b,c,d) from scalar products via Gram determinant.
# Ref: MertigBohmDenner1991, Eq. (2.21): ε(a,b,c,d)² = -det[G_{ij}]
function _evaluate_eps(e::Eps, sp::Dict{Tuple{Symbol,Symbol}, Float64})
    a, b, c, d = e.a::Momentum, e.b::Momentum, e.c::Momentum, e.d::Momentum
    vecs = (a.name, b.name, c.name, d.name)
    _sp(i, j) = sp[_sp_key(vecs[i], vecs[j])]
    G = [_sp(i,j) for i in 1:4, j in 1:4]
    det_G = det(G)
    abs_eps = sqrt(max(0.0, -det_G))
    _eps_sign(vecs) * abs_eps
end

# Sign convention for ε(names...) in standard 2→2 CM frame.
# Reference order: (:p1, :p2, :k1, :k2) has ε > 0 (sinθ > 0).
# Odd permutations flip sign. Unknown momenta → error.
function _eps_sign(names::NTuple{4,Symbol})
    ref = (:p1, :p2, :k1, :k2)
    perm = [findfirst(==(n), ref) for n in names]
    any(isnothing, perm) && error(
        "_eps_sign: unknown momentum name in $names; " *
        "only (:p1,:p2,:k1,:k2) supported for 2→2 sign convention")
    inv = sum(perm[i] > perm[j] for i in 1:4 for j in (i+1):4)
    iseven(inv) ? 1.0 : -1.0
end
