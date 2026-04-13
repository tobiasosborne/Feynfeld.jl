#  qg21 Step B: external-leg distribution (xc, xn enumeration).
#
#  Given a Partition + a TopoState whose vdeg layout is set, enumerate every
#  pair (xc, xn) that:
#
#    xc[1]  = 2 × #(external-to-external pairs)              [even, 0..n_ext]
#    xc[k]  = number of external legs going to internal      [k = 2..n1]
#             vertices in degree-class k
#    xn[i]  = number of external legs incident on vertex i
#             - For external vertices (i ≤ n_ext): always 0 (Step B never sets
#               externals' xn — qgraf treats them as singletons of degree 1).
#             - Within an internal degree-class, xn values are non-increasing
#               (canonicalisation: equivalent vertices ordered by leg count).
#             - Σ_{i ∈ class k} xn[i] = xc[k].
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12554-12658
#  Mapping (Fortran labels → Julia control flow):
#    label 61  →  outer xc-redistribution loop  (advance xc[2..n1])
#    label 43  →  initial xn fill within each class (greedy non-increasing)
#    label 76  →  xn-backtracking within innermost class
#    label 89  →  decrement xn[j] and repropagate
#    label 85  →  xc[1] += 2 (advance ext-to-ext count)
#    label 12  →  emit (callback) — handed back to caller (Step C)

"""
    _degree_class_bounds(vdeg, n) -> (xl::Vector{Int8}, xt::Vector{Int8})

Compute the degree-equivalence-class boundaries:
  xl[c]..xt[c] is the contiguous index range of vertices in class c
where a "class" is a maximal run of vertices sharing the same vdeg.
Includes external (vdeg=1) and internal classes alike.

Source: qgraf-4.0.6.f08:12541-12550.
"""
function _degree_class_bounds(vdeg::AbstractVector{Int8}, n::Integer)
    n = Int(n)
    xl = Int8[]
    xt = Int8[]
    n == 0 && return xl, xt
    push!(xl, Int8(1))
    @inbounds for i in 2:n
        if vdeg[i] != vdeg[i-1]
            push!(xt, Int8(i-1))
            push!(xl, Int8(i))
        end
    end
    push!(xt, Int8(n))
    return xl, xt
end

"""
    step_b_enumerate!(callback, state)

Enumerate every (xc, xn) external-leg distribution for the given Partition
encoded in `state`, calling `callback(state)` once per emission.  Mutates
`state.xc[1..n_int_classes+1]` and `state.xn[1..n]` in place; the callback
receives the same `state` (do not retain references between callbacks —
copy what you need).

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12554-12658.
Goto-style control flow mirrors the Fortran labels for line-by-line audit.
"""
function step_b_enumerate!(callback::F, state::TopoState) where {F}
    n      = Int(state.n)
    n_ext  = Int(state.n_ext)

    # qg21:12554-12557 — no externals: emit a single empty configuration.
    if n_ext == 0
        state.xc[1] = Int8(0)
        # Reset xn for safety (already zero from ctor, but clarify intent).
        @inbounds for i in 1:n
            state.xn[i] = Int8(0)
        end
        callback(state)
        return nothing
    end

    # Build degree-class boundaries.
    xl_all, xt_all = _degree_class_bounds(state.vdeg, n)

    # Identify internal classes (vdeg > 1).
    int_start = findfirst(c -> state.vdeg[xl_all[c]] > 1, eachindex(xl_all))
    if int_start === nothing
        # No internal vertices — all externals (process is trivial).
        # qgraf would still iterate xc(1)=0,2,4,...,n_ext, but with no internal
        # vertices the only valid case is xc(1)=n_ext (everything ext-to-ext).
        state.xc[1] = Int8(n_ext)
        @inbounds for i in 1:n
            state.xn[i] = Int8(0)
        end
        callback(state)
        return nothing
    end

    n_int_cls = length(xl_all) - int_start + 1

    # degr[j], bound[j] for j ∈ 1:n_int_cls (qgraf's degr/bound at index j+1).
    degr  = Int8[state.vdeg[xl_all[int_start + j - 1]] for j in 1:n_int_cls]
    bound = Int8[degr[j] * (xt_all[int_start + j - 1] - xl_all[int_start + j - 1] + Int8(1))
                 for j in 1:n_int_cls]
    kk    = sum(Int(b) for b in bound; init=0)

    # Per-internal-class vertex range from xl_all/xt_all
    xl(j) = Int(xl_all[int_start + j - 1])
    xt(j) = Int(xt_all[int_start + j - 1])

    # xc_int[j]  ↔  qgraf's xc(j+1)   (j ∈ 1:n_int_cls)
    xc_int = zeros(Int8, n_int_cls)
    xs     = zeros(Int8, n)
    j1     = 0          # current backtrack vertex index
    cur_cls = 0         # current backtrack class (1..n_int_cls)

    # qg21:12568 — initial xc(1)
    xc1 = Int8(max(0, n_ext - kk))

    @label outer_xc1
        # qg21:12569-12574 — greedy initial xc[2..n1] fill
        ii = n_ext - Int(xc1)
        @inbounds for j in n_int_cls:-1:1
            xc_int[j] = Int8(min(ii, Int(bound[j])))
            ii -= Int(xc_int[j])
        end
        # ii should be 0 if xc(1) is consistent.

    @label fill_xn
        # qg21:12576-12583 — greedy non-increasing xn fill within each class
        @inbounds for j in n_int_cls:-1:1
            ii = Int(xc_int[j])
            d  = Int(degr[j])
            for v in xl(j):xt(j)
                state.xn[v] = Int8(min(ii, d))
                ii -= Int(state.xn[v])
                xs[v] = Int8(ii)
            end
        end

    @label emit
        # Reflect xc into state.xc (state is the user-visible API).
        state.xc[1] = xc1
        @inbounds for j in 1:n_int_cls
            state.xc[j+1] = xc_int[j]
        end
        callback(state)

    # ── qg21:12585-12619 — xn backtrack (Fortran label 76, d01 loop) ───
    @label xn_backtrack
        # d01 loop: scan classes high → low looking for a vertex with xn>1.
        cur_cls = 0
        @inbounds for j in n_int_cls:-1:1
            ic_l = xl(j); ic_t = xt(j)
            for v in (ic_t - 1):-1:ic_l
                if state.xn[v] > Int8(1)
                    j1 = v
                    cur_cls = j
                    @goto decrement
                end
            end
        end
        # d01 exhausted — proceed to xc backtrack.
        @goto xc_backtrack

    @label decrement
        # qg21:12595-12609 — decrement xn[j1] and propagate.
        let j = cur_cls
            ic_l = xl(j); ic_t = xt(j)
            @inbounds begin
                state.xn[j1] -= Int8(1)
                xs[j1] += Int8(1)
                for v in (j1 + 1):ic_t
                    state.xn[v] = Int8(min(Int(state.xn[j1]), Int(xs[v - 1])))
                    xs[v] = xs[v - 1] - state.xn[v]
                end
                if xs[ic_t] > Int8(0)
                    # Inner backtrack: find smaller j1 in same class.
                    found = false
                    for v in (j1 - 1):-1:ic_l
                        if state.xn[v] > Int8(1)
                            j1 = v
                            found = true
                            break
                        end
                    end
                    if found
                        @goto decrement
                    else
                        # Cycle d01: try next class.
                        # Move to lower class and retry the d01 scan from there.
                        # We restart the d01 scan from class (cur_cls - 1) downward.
                        for jj in (cur_cls - 1):-1:1
                            ic_l2 = xl(jj); ic_t2 = xt(jj)
                            for v in (ic_t2 - 1):-1:ic_l2
                                if state.xn[v] > Int8(1)
                                    j1 = v
                                    cur_cls = jj
                                    @goto decrement
                                end
                            end
                        end
                        @goto xc_backtrack
                    end
                end
                # qg21:12610-12617 — reset lower classes to greedy default.
                for jj in (cur_cls + 1):n_int_cls
                    ic_l2 = xl(jj); ic_t2 = xt(jj)
                    ii2 = Int(xc_int[jj])
                    d2  = Int(degr[jj])
                    for v in ic_l2:ic_t2
                        state.xn[v] = Int8(min(ii2, d2))
                        ii2 -= Int(state.xn[v])
                        xs[v] = Int8(ii2)
                    end
                end
            end
        end
        @goto emit

    # ── qg21:12620-12651 — xc[2..n1] backtrack (Fortran label 45) ──────
    @label xc_backtrack
        # Find rightmost class j (Julia 1..n_int_cls; qgraf 3..n1) with xc>0.
        # Note: qgraf's "i1=n1..3 step -1" maps to Julia j=n_int_cls..2 step -1.
        let
            j1q = 0
            @inbounds for j in n_int_cls:-1:2
                if xc_int[j] > Int8(0)
                    j1q = j
                    break
                end
            end
            if j1q == 0
                @goto bump_xc1
            end
            # Find j2q ∈ (j1q-1)..1 (qgraf (j1q-1)..2 → 1-based here is j..1)
            # with xc_int[j] < bound[j].  qgraf's "i1=j1-1,2,-1" — its index 2
            # corresponds to our Julia index 1, so range is (j1q-1)..1.
            j2q = 0
            @inbounds for j in (j1q - 1):-1:1
                if xc_int[j] < bound[j]
                    j2q = j
                    break
                end
            end
            if j2q == 0
                @goto bump_xc1
            end
            # qg21:12636-12644 — increment and redistribute lower classes.
            @inbounds begin
                xc_int[j2q] += Int8(1)
                ii = -1
                for j in (j2q + 1):n_int_cls
                    ii += Int(xc_int[j])
                end
                for j in n_int_cls:-1:(j2q + 1)
                    xc_int[j] = Int8(min(ii, Int(bound[j])))
                    ii -= Int(xc_int[j])
                end
                ii == 0 || error("step_b: qg21_2 — invalid xc redistribution")
            end
        end
        @goto fill_xn

    @label bump_xc1
        # qg21:12652-12655 — xc(1) += 2.
        xc1 += Int8(2)
        if xc1 ≤ Int8(n_ext)
            @goto outer_xc1
        end
        return nothing
end
