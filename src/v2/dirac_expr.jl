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
# Keys on the FULL chain.elements (gammas + spinors). Two terms whose
# chains differ only in their spinor pair represent physically distinct
# fermion lines and must NOT merge.
# Bug history: feynfeld-ocpb (F002) — earlier this keyed on `gammas(chain)`,
# silently dropping spinors and collapsing distinct lines. Regression
# guard: test/v2/test_dirac_expr_simplify.jl.
function simplify(de::DiracExpr)
    groups = Dict{Vector{DiracElement}, AlgSum}()
    for (coeff, chain) in de.terms
        elems = chain.elements
        existing = get(groups, elems, AlgSum())
        groups[elems] = existing + coeff
    end
    terms = Tuple{AlgSum, DiracChain}[]
    for (elems, coeff) in groups
        iszero(coeff) && continue
        push!(terms, (coeff, DiracChain(elems)))
    end
    DiracExpr(terms)
end
