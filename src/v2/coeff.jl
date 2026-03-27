# DimPoly: polynomial in D for dimensional regularisation coefficients.
# Represents c₀ + c₁D + c₂D² + ... where cᵢ are Rational{Int}.
# Replaces SymDim (which could only handle aD+b and errored on D²).
# Future-proof: handles any power of D arising from traces, dim-reg, multi-loop.

struct DimPoly
    coeffs::Vector{Rational{Int}}  # coeffs[i] = coefficient of D^(i-1)
    # Invariant: no trailing zeros (trimmed), empty = zero
    function DimPoly(cs::Vector{Rational{Int}})
        # Trim trailing zeros
        n = length(cs)
        while n > 0 && iszero(cs[n])
            n -= 1
        end
        new(cs[1:n])
    end
end

DimPoly(cs::Vector{<:Integer}) = DimPoly(Rational{Int}.(cs))
DimPoly(cs::Vector{<:Number}) = DimPoly(Rational{Int}.(cs))

# Degree of the polynomial (-1 for zero polynomial)
degree(p::DimPoly) = length(p.coeffs) - 1

# Named constants
const DIM = DimPoly([0, 1])           # D
const DIM_MINUS_4 = DimPoly([-4, 1])  # D - 4

# Evaluate at a specific dimension value
evaluate_dim(p::DimPoly, d::Number=4) = isempty(p.coeffs) ? 0//1 : sum(c * d^(i-1) for (i,c) in enumerate(p.coeffs))
evaluate_dim(x::Number, d::Number=4) = x

# Pretty-print a rational as integer when possible
_cshow(c::Rational{Int}) = isinteger(c) ? string(Int(c)) : string(c)

function Base.show(io::IO, p::DimPoly)
    isempty(p.coeffs) && (print(io, "0"); return)
    # Collect non-zero terms, highest power first
    nonzero = [(i-1, c) for (i, c) in enumerate(p.coeffs) if !iszero(c)]
    reverse!(nonzero)  # highest power first
    first_term = true
    for (pow, c) in nonzero
        ac = abs(c)
        cs = _cshow(ac)
        if first_term
            neg = c < 0
            if pow == 0
                print(io, neg ? "-$cs" : cs)
            elseif isone(ac)
                print(io, neg ? "-D" : "D")
                pow > 1 && print(io, "^$pow")
            else
                print(io, neg ? "-$cs" : cs, "*D")
                pow > 1 && print(io, "^$pow")
            end
        else
            print(io, c > 0 ? " + " : " - ")
            if pow == 0
                print(io, cs)
            elseif isone(ac)
                print(io, "D")
                pow > 1 && print(io, "^$pow")
            else
                print(io, cs, "*D")
                pow > 1 && print(io, "^$pow")
            end
        end
        first_term = false
    end
end

# ---- Arithmetic: DimPoly <-> DimPoly ----
function Base.:+(a::DimPoly, b::DimPoly)
    n = max(length(a.coeffs), length(b.coeffs))
    cs = zeros(Rational{Int}, n)
    for (i, c) in enumerate(a.coeffs); cs[i] += c; end
    for (i, c) in enumerate(b.coeffs); cs[i] += c; end
    DimPoly(cs)
end

function Base.:-(a::DimPoly, b::DimPoly)
    n = max(length(a.coeffs), length(b.coeffs))
    cs = zeros(Rational{Int}, n)
    for (i, c) in enumerate(a.coeffs); cs[i] += c; end
    for (i, c) in enumerate(b.coeffs); cs[i] -= c; end
    DimPoly(cs)
end

Base.:-(p::DimPoly) = DimPoly(-p.coeffs)

function Base.:*(a::DimPoly, b::DimPoly)
    isempty(a.coeffs) && return a
    isempty(b.coeffs) && return b
    n = length(a.coeffs) + length(b.coeffs) - 1
    cs = zeros(Rational{Int}, n)
    for (i, ca) in enumerate(a.coeffs)
        for (j, cb) in enumerate(b.coeffs)
            cs[i + j - 1] += ca * cb
        end
    end
    DimPoly(cs)
end

Base.:(==)(a::DimPoly, b::DimPoly) = a.coeffs == b.coeffs
Base.hash(p::DimPoly, h::UInt) = hash(p.coeffs, hash(:DimPoly, h))
Base.iszero(p::DimPoly) = isempty(p.coeffs)
Base.isone(p::DimPoly) = length(p.coeffs) == 1 && isone(p.coeffs[1])
Base.zero(::Type{DimPoly}) = DimPoly(Rational{Int}[])
Base.one(::Type{DimPoly}) = DimPoly([1//1])

# ---- Arithmetic: DimPoly <-> Number ----
Base.:+(p::DimPoly, n::Number) = p + DimPoly([Rational{Int}(n)])
Base.:+(n::Number, p::DimPoly) = p + n
Base.:-(p::DimPoly, n::Number) = p - DimPoly([Rational{Int}(n)])
Base.:-(n::Number, p::DimPoly) = DimPoly([Rational{Int}(n)]) - p
Base.:*(p::DimPoly, n::Number) = DimPoly(p.coeffs .* Rational{Int}(n))
Base.:*(n::Number, p::DimPoly) = p * n

# ---- Coeff: the union type for coefficients ----
# 2-element union → Julia union-splitting is maximally efficient.
const Coeff = Union{Rational{Int}, DimPoly}

# Normalise: collapse DimPoly to Rational{Int} when it's a constant
function normalise_coeff(p::DimPoly)
    isempty(p.coeffs) && return 0//1
    length(p.coeffs) == 1 && return p.coeffs[1]
    p
end
normalise_coeff(x::Rational{Int}) = x
normalise_coeff(x::Integer) = Rational{Int}(x)
normalise_coeff(x::Number) = Rational{Int}(x)

# ---- mul_coeff / add_coeff: dispatch for Coeff operands ----
mul_coeff(a::Number, b::Number) = Rational{Int}(a) * Rational{Int}(b)
mul_coeff(a::DimPoly, b::DimPoly) = normalise_coeff(a * b)
mul_coeff(a::DimPoly, b::Number) = normalise_coeff(a * b)
mul_coeff(a::Number, b::DimPoly) = normalise_coeff(b * a)

add_coeff(a::Number, b::Number) = Rational{Int}(a) + Rational{Int}(b)
add_coeff(a::DimPoly, b::DimPoly) = normalise_coeff(a + b)
add_coeff(a::DimPoly, b::Number) = normalise_coeff(a + b)
add_coeff(a::Number, b::DimPoly) = normalise_coeff(b + a)
