# Scalar product lookup and K·p dot products for numerical evaluation.
# Used by TID (tid.jl) for evaluating loop integrals at specific kinematics.

# Symmetric lookup in Float64 scalar product dictionary.
function _splookup(sv::Dict{Tuple{Symbol,Symbol},Float64}, a::Symbol, b::Symbol)
    key = _sp_key(a, b)
    haskey(sv, key) && return sv[key]
    haskey(sv, (b, a)) && return sv[(b, a)]
    error("SP($a,$b) not in sp_vals")
end

# K · p_a where K is Momentum or MomentumSum, p_a is a Symbol
_Kdot(K::Momentum, pa::Symbol, sv) = _splookup(sv, K.name, pa)
function _Kdot(K::MomentumSum, pa::Symbol, sv)
    sum(Float64(c) * _splookup(sv, m.name, pa) for (c, m) in K.terms)
end

# K · K_ref where both are Momentum or MomentumSum
_Kdot(K::Momentum, K_ref::Momentum, sv) = _splookup(sv, K.name, K_ref.name)
function _Kdot(K::Momentum, K_ref::MomentumSum, sv)
    sum(Float64(c) * _splookup(sv, K.name, m.name) for (c, m) in K_ref.terms)
end
function _Kdot(K::MomentumSum, K_ref::Momentum, sv)
    sum(Float64(c) * _splookup(sv, m.name, K_ref.name) for (c, m) in K.terms)
end
function _Kdot(K::MomentumSum, K_ref::MomentumSum, sv)
    val = 0.0
    for (ci, mi) in K.terms, (cj, mj) in K_ref.terms
        val += Float64(ci) * Float64(cj) * _splookup(sv, mi.name, mj.name)
    end
    val
end
