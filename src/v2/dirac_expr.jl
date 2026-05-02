# DiracExpr: matrix-valued expressions in Dirac space.
# A DiracExpr is a sum of (AlgSum coefficient, DiracChain) pairs.
# This represents objects like the self-energy: Σ = A(p²)·I + B(p²)·p-slash
#
# Convention: DiracChain([]) = identity matrix (empty chain).

struct DiracExpr
    terms::Vector{Tuple{AlgSum, DiracChain}}
end

DiracExpr() = DiracExpr(Tuple{AlgSum, DiracChain}[])

Base.:(==)(a::DiracExpr, b::DiracExpr) = a.terms == b.terms
Base.hash(d::DiracExpr, h::UInt) = hash(d.terms, hash(:DiracExpr, h))

# Lift a scalar (AlgSum) to DiracExpr (scalar * identity)
DiracExpr(s::AlgSum) = DiracExpr([(s, DiracChain(DiracGamma[]))])
DiracExpr(n::Number) = DiracExpr(alg(n))

# Lift a DiracChain to DiracExpr (coefficient = 1)
DiracExpr(c::DiracChain) = DiracExpr([(alg(1), c)])

function Base.show(io::IO, de::DiracExpr)
    if isempty(de.terms)
        print(io, "0"); return
    end
    for (i, (coeff, chain)) in enumerate(de.terms)
        i > 1 && print(io, " + ")
        elems = chain.elements
        if isempty(elems)
            print(io, "(", coeff, ")·I")
        else
            print(io, "(", coeff, ")·[", chain, "]")
        end
    end
end

# ---- Arithmetic ----
function Base.:+(a::DiracExpr, b::DiracExpr)
    simplify(DiracExpr(vcat(a.terms, b.terms)))
end

function Base.:-(a::DiracExpr)
    DiracExpr([(-1 * c, chain) for (c, chain) in a.terms])
end

function Base.:-(a::DiracExpr, b::DiracExpr)
    a + (-b)
end

# Scalar * DiracExpr
function Base.:*(s::AlgSum, de::DiracExpr)
    DiracExpr([(s * c, chain) for (c, chain) in de.terms])
end
function Base.:*(de::DiracExpr, s::AlgSum)
    s * de
end
function Base.:*(n::Number, de::DiracExpr)
    alg(n) * de
end
function Base.:*(de::DiracExpr, n::Number)
    n * de
end
function Base.:*(d::DimPoly, de::DiracExpr)
    alg(d) * de
end
function Base.:*(de::DiracExpr, d::DimPoly)
    d * de
end

# DiracExpr * DiracExpr (matrix multiplication = chain concatenation)
function Base.:*(a::DiracExpr, b::DiracExpr)
    result = Tuple{AlgSum, DiracChain}[]
    for (ca, cha) in a.terms
        for (cb, chb) in b.terms
            push!(result, (ca * cb, DiracChain(vcat(cha.elements, chb.elements))))
        end
    end
    DiracExpr(result)
end

# ---- Take trace of DiracExpr → AlgSum ----
function dirac_trace(de::DiracExpr)
    result = AlgSum()
    for (coeff, chain) in de.terms
        mul_acc!(result, coeff, dirac_trace(chain))
    end
    result
end

# ---- Clean up: collect terms with identical chains ----
function simplify(de::DiracExpr)
    # Group by chain structure (using chain elements as key)
    groups = Dict{Vector{DiracGamma}, AlgSum}()
    for (coeff, chain) in de.terms
        gs = gammas(chain)
        existing = get(groups, gs, AlgSum())
        groups[gs] = existing + coeff
    end
    terms = Tuple{AlgSum, DiracChain}[]
    for (gs, coeff) in groups
        iszero(coeff) && continue
        push!(terms, (coeff, DiracChain(gs)))
    end
    DiracExpr(terms)
end
