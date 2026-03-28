# Passarino-Veltman loop integral functions.
# PaVe{N}: N-point scalar/tensor integral coefficient.
# Standalone type — not in AlgFactor union. Scalar-valued, no Lorentz indices.
# Named constructors: A0, B0, B1, C0, D0.

struct PaVe{N}
    indices::Vector{Int}        # tensor coefficient labels (sorted)
    invariants::Vector{Float64} # N*(N-1)/2 kinematic invariants
    masses::Vector{Float64}     # N squared masses

    function PaVe{N}(indices::Vector{Int},
                     invariants::Vector{Float64},
                     masses::Vector{Float64}) where {N}
        N >= 1 || throw(ArgumentError("PaVe{N}: N must be >= 1, got $N"))
        n_inv = N * (N - 1) ÷ 2
        length(invariants) == n_inv ||
            throw(ArgumentError("PaVe{$N}: need $n_inv invariants, got $(length(invariants))"))
        length(masses) == N ||
            throw(ArgumentError("PaVe{$N}: need $N masses, got $(length(masses))"))
        new{N}(sort(indices), invariants, masses)
    end
end

# ---- Named constructors ----

A0(m2::Real) = PaVe{1}(Int[], Float64[], [Float64(m2)])

B0(p2::Real, m02::Real, m12::Real) =
    PaVe{2}(Int[], [Float64(p2)], [Float64(m02), Float64(m12)])

B1(p2::Real, m02::Real, m12::Real) =
    PaVe{2}([1], [Float64(p2)], [Float64(m02), Float64(m12)])

B00(p2::Real, m02::Real, m12::Real) =
    PaVe{2}([0, 0], [Float64(p2)], [Float64(m02), Float64(m12)])

B11(p2::Real, m02::Real, m12::Real) =
    PaVe{2}([1, 1], [Float64(p2)], [Float64(m02), Float64(m12)])

C0(p10::Real, p12::Real, p20::Real, m02::Real, m12::Real, m22::Real) =
    PaVe{3}(Int[], [Float64(p10), Float64(p12), Float64(p20)],
             [Float64(m02), Float64(m12), Float64(m22)])

C1(p10::Real, p12::Real, p20::Real, m02::Real, m12::Real, m22::Real) =
    PaVe{3}([1], [Float64(p10), Float64(p12), Float64(p20)],
             [Float64(m02), Float64(m12), Float64(m22)])

C2(p10::Real, p12::Real, p20::Real, m02::Real, m12::Real, m22::Real) =
    PaVe{3}([2], [Float64(p10), Float64(p12), Float64(p20)],
             [Float64(m02), Float64(m12), Float64(m22)])

D0(p10::Real, p12::Real, p23::Real, p30::Real, p20::Real, p13::Real,
   m02::Real, m12::Real, m22::Real, m32::Real) =
    PaVe{4}(Int[], [Float64(p10), Float64(p12), Float64(p23),
                     Float64(p30), Float64(p20), Float64(p13)],
             [Float64(m02), Float64(m12), Float64(m22), Float64(m32)])

# ---- Standard Julia interface ----

function Base.:(==)(a::PaVe{N}, b::PaVe{N}) where {N}
    a.indices == b.indices && a.invariants == b.invariants && a.masses == b.masses
end
Base.:(==)(::PaVe{N}, ::PaVe{M}) where {N, M} = false

function Base.hash(p::PaVe{N}, h::UInt) where {N}
    hash(p.masses, hash(p.invariants, hash(p.indices, hash(N, hash(:PaVe, h)))))
end

const _PAVE_NAMES = Dict(1 => "A", 2 => "B", 3 => "C", 4 => "D")

function Base.show(io::IO, p::PaVe{N}) where {N}
    name = get(_PAVE_NAMES, N, "PaVe{$N}")
    idx_str = join(p.indices, "")
    args = join(vcat(p.invariants, p.masses), ", ")
    print(io, name, idx_str, "(", args, ")")
end
