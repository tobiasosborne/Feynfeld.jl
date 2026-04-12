# Strategy C: hybrid recursive-descent port of qgraf's qg21 topology generator.
#
# Ref: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08 (Paulo Nogueira)
# Ref: refs/papers/Nogueira1993_JCompPhys105_279.pdf §2 "The method"
# Ref: refs/qgraf/ALGORITHM.md (cleanroom spec)
#
# Beads: feynfeld-ney (master), feynfeld-5hr (this phase).
#
# No StaticArrays dependency. Scratch arrays are preallocated to MAX_V size
# on TopoState construction; the recursive inner loop reads and writes in
# place so the hot path allocates zero per backtrack frame.

"""
    MAX_V

Maximum total vertex count (external + internal) supported by TopoState
scratch arrays. qgraf's default is 20ish; Feynfeld targets practical QFT
processes at ≤ 4 loops where n ≤ ~18.
"""
const MAX_V = 24

"""
    Partition(n_ext, counts, mrho, nloop)

A vertex-degree partition ρ, as in qgraf's qpg11.

- `n_ext`: ρ(-1), number of external legs (each of degree 1).
- `counts`: counts[i] = number of internal vertices of degree `mrho + i - 1`
  for i ∈ 1:length(counts). So `counts[1]` is ρ(mrho), ..., `counts[end]` is ρ(nrho).
- `mrho`: smallest internal degree present in the partition.
- `nloop`: prescribed loop count.

Constraints (Euler + degree balance, ALGORITHM.md §2.1):
  Σ_k k·ρ(k) + ρ(-1) = 2·P            (total half-edges)
  Σ_k ρ(k) + ρ(-1)   = V               (total vertices)
  L = P − V + 1                         (first Betti number)

Ref: qgraf-4.0.6.f08:1743-1991 (qpg11), ALGORITHM.md §2.
"""
struct Partition
    n_ext::Int8
    counts::Vector{Int8}   # ρ(mrho..nrho)
    mrho::Int8
    nloop::Int8
end

"Count ρ(k) of internal vertices of degree `k` (0 if out of range)."
function rho_k(p::Partition, k::Integer)::Int8
    idx = Int(k) - Int(p.mrho) + 1
    (idx < 1 || idx > length(p.counts)) && return Int8(0)
    p.counts[idx]
end

"Total number of internal vertices: Σ_k ρ(k)."
n_internal(p::Partition) = Int(sum(p.counts; init=Int8(0)))

"Total vertices V = ρ(-1) + Σ_k ρ(k)."
n_vertices(p::Partition) = Int(p.n_ext) + n_internal(p)

"Total internal edges P from Σ_k k·ρ(k) + ρ(-1) = 2·P."
function n_edges(p::Partition)::Int
    s = Int(p.n_ext)
    for (i, c) in enumerate(p.counts)
        s += (Int(p.mrho) + i - 1) * Int(c)
    end
    s ÷ 2
end

"""
    EquivClass(first, last, degree)

A contiguous range of internal vertices sharing the same (degree, xn, xg_diag).
Vertices are numbered from 1, with externals at 1..n_ext, then internals in
ascending (degree, xn, self-loop-count) order. qgraf's `uset`/`xset` arrays
encode exactly this partition.

Ref: qg21 labels 19/28 (xset/uset construction, qgraf-4.0.6.f08:12819-12840).
"""
struct EquivClass
    members::Vector{Int8}   # vertex indices, sorted ascending; not necessarily contiguous
    degree::Int8
end

# Convenience: build from a contiguous range (used in qg21 Step A/B where
# vertices are guaranteed sorted by invariants).
function EquivClass(first::Int8, last::Int8, degree::Int8)
    first <= last + Int8(1) || error("EquivClass: first=$first > last+1=$(last+1)")
    members = Int8[Int8(i) for i in Int(first):Int(last)]
    EquivClass(members, degree)
end

Base.length(c::EquivClass) = length(c.members)
Base.iterate(c::EquivClass, state::Int=1) =
    state > length(c.members) ? nothing : (c.members[state], state + 1)
Base.first(c::EquivClass) = c.members[1]
Base.last(c::EquivClass) = c.members[end]

"""
    FilterSet(; kwargs...)

Topological filters mirroring qgraf's dflag option space.
Ref: qgraf-4.0.6.f08:3690 (qumpi), 3777 (qumvi), 13669 (qgsig), 18830 (qcyc).
Ref: ALGORITHM.md §3.7.
"""
Base.@kwdef struct FilterSet
    onepi::Bool       = false
    nobridge::Bool    = false
    nosbridge::Bool   = false
    notadpole::Bool   = false
    onshell::Bool     = false
    nosnail::Bool     = false
    onevi::Bool       = false
    onshellx::Bool    = false
    noselfloop::Bool  = false
    nodiloop::Bool    = false
    noparallel::Bool  = false
    cycli::Bool       = false
    nosigma::Bool     = false
    bipart::Bool      = false
end

"True iff no filter is enabled — the default construct `FilterSet()`."
no_filters(f::FilterSet) = !(f.onepi | f.nobridge | f.nosbridge | f.notadpole |
                             f.onshell | f.nosnail | f.onevi | f.onshellx |
                             f.noselfloop | f.nodiloop | f.noparallel |
                             f.cycli | f.nosigma | f.bipart)

"""
    TopoState(partition)

Mutable scratch state for qg21's row-by-row topology construction.
All arrays preallocated to MAX_V size; hot path allocates nothing.
Phase 1 subset of qg21 SAVEd locals (qgraf-4.0.6.f08:12426-12446):
xg/xn/xc/ds/lps/dta/xset/classes/ngsym + perm/aa/queue scratch.
uset, str, xl, xt, xp, a1, bb, p1s, head/tail/intree, emul/nemul are
deferred to later phases as they come into use.
"""
mutable struct TopoState
    n::Int8
    n_ext::Int8
    rhop1::Int8
    nloop::Int8
    vdeg::Vector{Int8}
    xg::Matrix{Int8}
    xn::Vector{Int8}
    xc::Vector{Int8}
    ds::Matrix{Int8}
    lps::Vector{Int8}
    dta::Vector{Int8}
    xset::Vector{Int8}
    classes::Vector{EquivClass}
    ngsym::Int32
    perm_buf::Vector{Int8}
    aa_buf::Vector{Int8}
    queue_buf::Vector{Int8}
end

function TopoState(p::Partition)
    n = n_vertices(p)
    n ≤ MAX_V || error("TopoState: n=$n exceeds MAX_V=$MAX_V")
    vdeg = zeros(Int8, MAX_V)
    # Externals at 1..n_ext get degree 1.
    for i in 1:Int(p.n_ext)
        vdeg[i] = Int8(1)
    end
    # Internals at n_ext+1..n get degrees in ascending order of `counts`.
    pos = Int(p.n_ext) + 1
    for (i, c) in enumerate(p.counts)
        d = Int(p.mrho) + i - 1
        for _ in 1:Int(c)
            vdeg[pos] = Int8(d)
            pos += 1
        end
    end
    TopoState(
        Int8(n),
        Int8(p.n_ext),
        Int8(p.n_ext + 1),
        p.nloop,
        vdeg,
        zeros(Int8, MAX_V, MAX_V),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V, MAX_V),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V),
        EquivClass[],
        Int32(0),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V),
        zeros(Int8, MAX_V),
    )
end
