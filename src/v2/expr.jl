# AlgExpr: expression tree for Lorentz algebra results.
# Design: Dict-based AlgSum (SymbolicUtils Add pattern) for O(1) like-term collection.
# Coefficients are Coeff = Union{Rational{Int}, DimPoly} — never Any, never Expr.

# ---- Factor types ----
struct Eps
    a::PairArg
    b::PairArg
    c::PairArg
    d::PairArg
end
Base.hash(e::Eps, h::UInt) = hash(e.a, hash(e.b, hash(e.c, hash(e.d, hash(:Eps, h)))))
Base.:(==)(a::Eps, b::Eps) = a.a == b.a && a.b == b.b && a.c == b.c && a.d == b.d

# The union of things that can appear as factors in an algebraic term.
# 6-element union: Lorentz (Pair, Eps) + Colour (SUNDelta, FundDelta, SUNF, SUND)
const AlgFactor = Union{Pair, Eps, SUNDelta, FundDelta, SUNF, SUND}

# ---- Structural ordering for factors (replaces repr-based sorting) ----
_factor_type_tag(::Pair) = 1
_factor_type_tag(::Eps) = 2
_factor_type_tag(::SUNDelta) = 3
_factor_type_tag(::FundDelta) = 4
_factor_type_tag(::SUNF) = 5
_factor_type_tag(::SUND) = 6

function Base.isless(a::Eps, b::Eps)
    a.a != b.a && return isless(a.a, b.a)
    a.b != b.b && return isless(a.b, b.b)
    a.c != b.c && return isless(a.c, b.c)
    return isless(a.d, b.d)
end

function _factor_isless(a::AlgFactor, b::AlgFactor)
    ta, tb = _factor_type_tag(a), _factor_type_tag(b)
    ta != tb && return ta < tb
    isless(a, b)
end

# ---- FactorKey: canonical sorted factor list for Dict lookup ----
struct FactorKey
    factors::Vector{AlgFactor}
    function FactorKey(fs::Vector{<:AlgFactor})
        sorted = sort(collect(AlgFactor, fs); lt=_factor_isless)
        new(sorted)
    end
end
FactorKey() = FactorKey(AlgFactor[])

Base.:(==)(a::FactorKey, b::FactorKey) = a.factors == b.factors
Base.hash(fk::FactorKey, h::UInt) = hash(fk.factors, hash(:FactorKey, h))

# ---- AlgSum: Dict{FactorKey, Coeff} ----
# Invariant: no zero-coefficient entries. Empty dict = zero.
# Values are Coeff (Rational{Int} or DimPoly) — typed, not Any.
struct AlgSum
    terms::Dict{FactorKey, Coeff}
end

AlgSum() = AlgSum(Dict{FactorKey, Coeff}())

function Base.show(io::IO, s::AlgSum)
    if isempty(s.terms)
        print(io, "0"); return
    end
    first_term = true
    for (fk, c) in s.terms
        first_term || print(io, " + ")
        first_term = false
        if isempty(fk.factors)
            print(io, c)
        elseif c isa Number && isone(c)
            join(io, fk.factors, "*")
        else
            print(io, c, "*")
            join(io, fk.factors, "*")
        end
    end
end

Base.iszero(s::AlgSum) = isempty(s.terms)

# ---- Constructors ----
_d() = Dict{FactorKey, Coeff}
alg(f::AlgFactor) = AlgSum(_d()(FactorKey([f]) => 1//1))
alg(c::Number) = iszero(c) ? AlgSum() : AlgSum(_d()(FactorKey() => Rational{Int}(c)))
alg(c::DimPoly) = iszero(c) ? AlgSum() : AlgSum(_d()(FactorKey() => c))
alg(s::AlgSum) = s

# ---- Helper: is a Coeff zero? ----
_coeff_iszero(c::Rational{Int}) = iszero(c)
_coeff_iszero(c::DimPoly) = iszero(c)
_coeff_iszero(c::Number) = iszero(c)

# ---- Arithmetic on AlgSum ----
function Base.:+(a::AlgSum, b::AlgSum)
    result = _d()(a.terms)
    for (fk, c) in b.terms
        existing = get(result, fk, 0//1)
        new_c = normalise_coeff(add_coeff(existing, c))
        if _coeff_iszero(new_c)
            delete!(result, fk)
        else
            result[fk] = new_c
        end
    end
    AlgSum(result)
end

Base.:-(a::AlgSum, b::AlgSum) = a + (-1 * b)
Base.:-(a::AlgSum) = -1 * a

function Base.:*(a::AlgSum, b::AlgSum)
    result = _d()()
    for (fk_a, c_a) in a.terms
        for (fk_b, c_b) in b.terms
            merged = FactorKey(vcat(fk_a.factors, fk_b.factors))
            c = normalise_coeff(mul_coeff(c_a, c_b))
            _coeff_iszero(c) && continue
            existing = get(result, merged, 0//1)
            new_c = normalise_coeff(add_coeff(existing, c))
            if _coeff_iszero(new_c)
                delete!(result, merged)
            else
                result[merged] = new_c
            end
        end
    end
    AlgSum(result)
end

# Scalar * AlgSum
function Base.:*(c::Number, s::AlgSum)
    iszero(c) && return AlgSum()
    isone(c) && return s
    AlgSum(_d()(fk => normalise_coeff(mul_coeff(v, c)) for (fk, v) in s.terms))
end
Base.:*(s::AlgSum, c::Number) = c * s
function Base.:*(c::DimPoly, s::AlgSum)
    iszero(c) && return AlgSum()
    AlgSum(_d()(fk => normalise_coeff(mul_coeff(v, c)) for (fk, v) in s.terms))
end
Base.:*(s::AlgSum, c::DimPoly) = c * s

# Number + AlgSum
Base.:+(c::Number, s::AlgSum) = alg(c) + s
Base.:+(s::AlgSum, c::Number) = s + alg(c)
Base.:-(c::Number, s::AlgSum) = alg(c) - s
Base.:-(s::AlgSum, c::Number) = s - alg(c)
