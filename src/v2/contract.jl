# Contraction engine: Einstein summation via Lorentz index matching.
# Returns AlgSum always — no mixed return types.

function contract(s::AlgSum; ctx::SPContext=CURRENT_SP[])
    result = AlgSum()
    for (fk, c) in s.terms
        contracted = _contract_factors(fk.factors, c, ctx)
        result = result + contracted
    end
    result
end

function _contract_factors(factors::Vector{AlgFactor}, coeff, ctx::SPContext)
    factors = copy(factors)
    coeff_acc = coeff

    # Worklist: repeatedly scan for contractable index pairs
    changed = true
    while changed
        changed = false
        indices = _index_inventory(factors)

        for (idx, positions) in indices
            length(positions) == 2 || continue

            i, j = positions
            if i == j
                # Same factor has both indices (trace): g^mu_mu = D
                result = _do_self_contraction(factors[i], idx)
            else
                fi, fj = factors[i], factors[j]
                result = _do_contraction(fi, fj, idx, ctx)
            end

            if result !== nothing
                new_coeff, new_factors = result
                # Remove old factors (largest index first)
                if i == j
                    deleteat!(factors, i)
                else
                    hi, lo = max(i,j), min(i,j)
                    deleteat!(factors, hi)
                    deleteat!(factors, lo)
                end
                append!(factors, new_factors)
                coeff_acc = mul_coeff(coeff_acc, new_coeff)
                coeff_acc = normalise_coeff(coeff_acc)
                changed = true
                break
            end
        end
    end

    alg_from_factors(factors, coeff_acc)
end

# Build an AlgSum from factors and coefficient
function alg_from_factors(factors::Vector{AlgFactor}, coeff)
    coeff = normalise_coeff(coeff)
    _coeff_iszero(coeff) && return AlgSum()
    AlgSum(Dict{FactorKey, Coeff}(FactorKey(factors) => coeff))
end

# Index inventory: which indices appear where?
function _index_inventory(factors::Vector{AlgFactor})
    inventory = Dict{LorentzIndex, Vector{Int}}()
    for (i, f) in enumerate(factors)
        for idx in _indices(f)
            push!(get!(inventory, idx, Int[]), i)
        end
    end
    inventory
end

# Extract Lorentz indices from a factor (dispatch on Pair type parameters!)
_indices(p::Pair{LorentzIndex, LorentzIndex}) = [p.a, p.b]
_indices(p::Pair{LorentzIndex, Momentum}) = [p.a]
_indices(p::Pair{Momentum, LorentzIndex}) = [p.b]
_indices(p::Pair{Momentum, Momentum}) = LorentzIndex[]
_indices(::Eps) = LorentzIndex[]  # TODO: for Eps
_indices(::AlgFactor) = LorentzIndex[]

# ---- Self-contraction: same factor has both occurrences of idx ----
function _do_self_contraction(f::MetricTensor, idx::LorentzIndex)
    f.a == idx && f.b == idx || return nothing
    (dim_trace(f.a.dim), AlgFactor[])
end
_do_self_contraction(::AlgFactor, ::LorentzIndex) = nothing

# ---- Contraction between two distinct factors ----

# Metric-Metric: g^{mu nu} g_{mu rho} → g^{nu rho}
function _do_contraction(fi::MetricTensor, fj::MetricTensor, idx::LorentzIndex, ::SPContext)
    surv_i = _surviving(fi, idx)
    surv_j = _surviving(fj, idx)
    surv_i !== nothing && surv_j !== nothing || return nothing
    dc = dim_contract(surv_i.dim, surv_j.dim)
    dc === nothing && return (0, AlgFactor[])
    new_p = pair(surv_i, surv_j)
    new_p isa Number && return (new_p, AlgFactor[])
    (1, AlgFactor[new_p])
end

# Metric-FourVector: g^{mu nu} p^mu → p^nu
function _do_contraction(fi::MetricTensor, fj::Pair{LorentzIndex, Momentum}, idx::LorentzIndex, ::SPContext)
    surv = _surviving(fi, idx)
    surv !== nothing || return nothing
    fj.a == idx || return nothing
    new_p = pair(surv, fj.b)
    new_p isa Number && return (new_p, AlgFactor[])
    (1, AlgFactor[new_p])
end

function _do_contraction(fi::Pair{LorentzIndex, Momentum}, fj::MetricTensor, idx::LorentzIndex, ctx::SPContext)
    _do_contraction(fj, fi, idx, ctx)
end

# FourVector-FourVector: p^mu q^mu → p.q
function _do_contraction(fi::Pair{LorentzIndex, Momentum}, fj::Pair{LorentzIndex, Momentum}, idx::LorentzIndex, ctx::SPContext)
    fi.a == idx && fj.a == idx || return nothing
    sp = pair(fi.b, fj.b)
    sp isa Number && return (sp, AlgFactor[])
    val = _try_sp(sp, ctx)
    val !== nothing && return (val, AlgFactor[])
    (1, AlgFactor[sp])
end

# Fallback
_do_contraction(::AlgFactor, ::AlgFactor, ::LorentzIndex, ::SPContext) = nothing

# ---- Index substitution (for polarization sums) ----
# Replace every occurrence of `old_idx` with `new_idx` in an AlgSum.

function substitute_index(s::AlgSum, old_idx::LorentzIndex, new_idx::LorentzIndex)
    result = AlgSum()
    for (fk, c) in s.terms
        new_factors = AlgFactor[_subst_factor(f, old_idx, new_idx) for f in fk.factors]
        result = result + alg_from_factors(new_factors, c)
    end
    result
end

_subst_li(li::LorentzIndex, old::LorentzIndex, new::LorentzIndex) =
    li == old ? new : li
_subst_mom(m::Momentum, ::LorentzIndex, ::LorentzIndex) = m
_subst_mom(m::MomentumSum, ::LorentzIndex, ::LorentzIndex) = m
_subst_pairarg(pa::LorentzIndex, old, new) = _subst_li(pa, old, new)
_subst_pairarg(pa::Momentum, old, new) = pa
_subst_pairarg(pa::MomentumSum, old, new) = pa

function _subst_factor(p::Pair, old::LorentzIndex, new::LorentzIndex)
    pair(_subst_pairarg(p.a, old, new), _subst_pairarg(p.b, old, new))
end
_subst_factor(f::AlgFactor, ::LorentzIndex, ::LorentzIndex) = f  # non-Pair: no indices to substitute

# Surviving index: which index of a MetricTensor is NOT the contracted one?
function _surviving(p::MetricTensor, idx::LorentzIndex)
    p.a == idx && return p.b
    p.b == idx && return p.a
    nothing
end
