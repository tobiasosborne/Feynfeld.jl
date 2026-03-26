# Feynfeld.jl — AlgTerm/AlgSum: algebraic expression tree
#
# Canonical sum-of-products representation for Lorentz algebra results.
# AlgTerm = coeff × ∏ factors (Pair, Eps — commutative)
# AlgSum  = Σ AlgTerms
#
# All algebra operations (contract, dirac_trace, expand_scalar_product)
# compose through AlgSum. Arithmetic: +, *, scalar multiply.
#
# This does NOT replace FTerm/Amplitude (which represent Feynman amplitudes).

export AlgTerm, AlgSum, alg, alg_zero, alg_scalar, is_scalar

# ── Factor type ────────────────────────────────────────────────────

"""Commutative factors that carry Lorentz indices."""
const AlgFactor = Union{Feynfeld.Pair, Eps}

# ── AlgTerm ────────────────────────────────────────────────────────

"""
    AlgTerm(coeff, factors)

A single monomial: `coeff × factor₁ × factor₂ × ⋯`
where each factor is a Pair or Eps (commutative product).
Coefficient can be Number, Symbol, or Expr.
"""
struct AlgTerm
    coeff::Any
    factors::Vector{AlgFactor}
end

AlgTerm(c) = AlgTerm(c, AlgFactor[])
AlgTerm(c, f::AlgFactor) = AlgTerm(c, AlgFactor[f])

# ── AlgSum ─────────────────────────────────────────────────────────

"""
    AlgSum(terms)

A sum of AlgTerms: `term₁ + term₂ + ⋯`
Canonical form for all algebraic computation results.
"""
struct AlgSum <: FeynExpr
    terms::Vector{AlgTerm}
end

AlgSum() = AlgSum(AlgTerm[])
AlgSum(t::AlgTerm) = AlgSum([t])

# ── Constructors / lifters ─────────────────────────────────────────

"""Lift a scalar, Pair, or Eps into an AlgSum."""
alg(x::Number) = AlgSum(AlgTerm[AlgTerm(x)])
alg(x::Symbol) = AlgSum(AlgTerm[AlgTerm(x)])
alg(x::Expr) = AlgSum(AlgTerm[AlgTerm(x)])
alg(p::Feynfeld.Pair) = AlgSum(AlgTerm[AlgTerm(1, AlgFactor[p])])
alg(e::Eps) = AlgSum(AlgTerm[AlgTerm(1, AlgFactor[e])])
alg(s::AlgSum) = s

"""Zero expression."""
alg_zero() = AlgSum(AlgTerm[])

"""Scalar expression (no Pair/Eps factors)."""
alg_scalar(c) = AlgSum(AlgTerm[AlgTerm(c)])

"""Check if an AlgSum is purely scalar (no Lorentz-indexed factors)."""
function is_scalar(s::AlgSum)
    all(t -> isempty(t.factors), s.terms)
end

# ── Factor ordering (structural, not repr-based) ──────────────────

_factor_key(p::Feynfeld.Pair) = (0, repr(p.a), repr(p.b))
_factor_key(e::Eps) = (1, repr(e.args))

function _sorted_factors(t::AlgTerm)
    sort(t.factors; by=_factor_key)
end

# ── Equality ───────────────────────────────────────────────────────

function Base.:(==)(a::AlgTerm, b::AlgTerm)
    a.coeff == b.coeff && _sorted_factors(a) == _sorted_factors(b)
end

Base.hash(t::AlgTerm, h::UInt) = hash(_sorted_factors(t), hash(t.coeff, hash(:AlgTerm, h)))

function Base.:(==)(a::AlgSum, b::AlgSum)
    ca, cb = _collect_terms(a), _collect_terms(b)
    length(ca.terms) != length(cb.terms) && return false
    Set(ca.terms) == Set(cb.terms)
end

function Base.hash(s::AlgSum, h::UInt)
    c = _collect_terms(s)
    hash(Set(c.terms), hash(:AlgSum, h))
end

# ── Arithmetic: addition ───────────────────────────────────────────

function Base.:+(a::AlgSum, b::AlgSum)
    _collect_terms(AlgSum(vcat(a.terms, b.terms)))
end

Base.:+(a::AlgSum, b::Number) = a + alg(b)
Base.:+(a::Number, b::AlgSum) = alg(a) + b
Base.:+(a::AlgSum, b::Feynfeld.Pair) = a + alg(b)
Base.:+(a::Feynfeld.Pair, b::AlgSum) = alg(a) + b

Base.:-(a::AlgSum) = AlgSum([AlgTerm(_mul_coeff(-1, t.coeff), copy(t.factors)) for t in a.terms])
Base.:-(a::AlgSum, b::AlgSum) = a + (-b)

# ── Arithmetic: multiplication ─────────────────────────────────────

function Base.:*(a::AlgSum, b::AlgSum)
    terms = AlgTerm[]
    for ta in a.terms, tb in b.terms
        c = _mul_coeff(ta.coeff, tb.coeff)
        c == 0 && continue
        fs = vcat(ta.factors, tb.factors)
        push!(terms, AlgTerm(c, fs))
    end
    _collect_terms(AlgSum(terms))
end

Base.:*(a::AlgSum, b::Number) = _scale(a, b)
Base.:*(a::Number, b::AlgSum) = _scale(b, a)
Base.:*(a::AlgSum, b::Symbol) = _scale(a, b)
Base.:*(a::Symbol, b::AlgSum) = _scale(b, a)
Base.:*(a::AlgSum, b::Feynfeld.Pair) = a * alg(b)
Base.:*(a::Feynfeld.Pair, b::AlgSum) = alg(a) * b

function _scale(s::AlgSum, c)
    c == 0 && return alg_zero()
    c == 1 && return s
    AlgSum([AlgTerm(_mul_coeff(t.coeff, c), copy(t.factors)) for t in s.terms])
end

# ── Collecting like terms ──────────────────────────────────────────

"""Collect terms with identical factor sets, summing their coefficients."""
function _collect_terms(s::AlgSum)
    isempty(s.terms) && return s
    groups = Dict{Vector{AlgFactor},Any}()
    for t in s.terms
        t.coeff == 0 && continue
        sf = _sorted_factors(t)
        if haskey(groups, sf)
            groups[sf] = _add_coeff(groups[sf], t.coeff)
        else
            groups[sf] = t.coeff
        end
    end
    terms = AlgTerm[AlgTerm(c, fs) for (fs, c) in groups if c != 0]
    AlgSum(terms)
end

# ── Coefficient arithmetic ─────────────────────────────────────────

function _add_coeff(a, b)
    a == 0 && return b
    b == 0 && return a
    (a isa Number && b isa Number) && return a + b
    :($a + $b)
end
