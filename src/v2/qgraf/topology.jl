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

"""
    qg21_enumerate!(callback, state)

Drive the full qg21 topology enumeration: Step B (xc/xn distribution) feeds
each emission into Step C (xg topology generation).  Calls `callback(state)`
once per emitted canonical topology.

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12426-13150.
"""
function qg21_enumerate!(callback::F, state::TopoState) where {F}
    step_b_enumerate!(state) do s
        step_c_enumerate!(callback, s)
    end
    return nothing
end

"""
    step_c_enumerate!(callback, state)

Step C: given an (xc, xn) configuration in `state`, enumerate all
non-isomorphic topologies xg respecting vdeg, xn and xc.  Calls
`callback(state)` once per emitted topology (canonical adjacency matrix
in `state.xg[1:n, 1:n]`).

Source: qgraf-4.0.6.f08:12659-13150.

Goto-style control flow (Fortran labels mapped 1:1 to Julia @label/@goto):

  setup       (12695-12764)  ← compute str, dta, place ext edges
  label 70    (12765-12781)  ← outer dsum loop, xg-diagonal init
  label 19    (12815-12841)  ← ds[] init, xset build, lps[lin]=loop-dsum+1
  label 10    (12842-12853)  ← greedy row fill at current lin
  label 28    (12897-12958)  ← row-completion check, advance to next row or emit
  label 200   (12959+)       ← emit topology
  label 17    (12862-12876)  ← row backtrack (decrement and re-fill)
  label 23    (12877-12890)  ← column-increment within row
  label 38    (12891-12896)  ← cross-col canonical-violation re-entry to 23
  label 15    (12857-12861)  ← row-decrement entry
  label 32    (12782-12798)  ← xg-diagonal backtrack (find rightmost xg(i,i)>0)
  label 33    (12790-12798)  ← find smaller-index slot to bump
  label 65    (12799-12814)  ← redistribute lower xg(i,i) after bump
"""
function step_c_enumerate!(callback::F, state::TopoState) where {F}
    n      = Int(state.n)
    n_ext  = Int(state.n_ext)
    rhop1  = Int(state.rhop1)
    nloop  = Int(state.nloop)

    # qg21:12660-12664 — reject xc(1)>0 when rho(-1)<n.
    if n_ext < n && state.xc[1] > Int8(0)
        return nothing
    end

    # ── Pre-compute str[], dta[] (default filters only — Phase 14 will
    #     extend for nob/nosb/nopa/etc.) ────────────────────────────────
    str_arr = zeros(Int8, MAX_V)
    dta_arr = zeros(Int8, MAX_V)
    @inbounds for i in rhop1:n
        # qg21:12530-12540 — str[i] = max edges from row i
        if state.vdeg[i] != state.vdeg[n]
            str_arr[i] = Int8(min(Int(state.vdeg[i]), nloop + 1))
        elseif n > 2
            str_arr[i] = Int8(min(Int(state.vdeg[i]) - 1, nloop + 1))
        else
            str_arr[i] = Int8(min(Int(state.vdeg[i]), nloop + 1))
        end
        # qg21:12722-12740 — dta[i] = max self-loops × 2 (in half-edges)
        ii_dta = (n > rhop1) ? 1 : 0
        jj = max(0, min(Int(state.vdeg[i]) - Int(state.xn[i]) - ii_dta, 2 * nloop)) ÷ 2
        dta_arr[i] = Int8(2 * jj)
    end

    # ── qg21:12743-12764 — clear xg, place ext-to-ext + ext-to-internal ─
    @inbounds for i in 1:n
        for j in 1:n
            state.xg[i, j] = Int8(0)
        end
    end
    half_pairs = Int(state.xc[1]) ÷ 2
    @inbounds for i in 1:half_pairs
        state.xg[2i - 1, 2i] = Int8(1)
    end
    jcur = Int(state.xc[1])
    @inbounds for i in rhop1:n
        if state.xn[i] > Int8(0)
            for _ in 1:Int(state.xn[i])
                jcur += 1
                state.xg[jcur, i] = Int8(1)
            end
        end
    end

    # ── State machine variables ─────────────────────────────────────────
    limin = rhop1
    limax = max(limin, n - 1)
    lin   = limin
    dsum  = -1
    aux   = 0
    col   = 0
    msum  = 0
    j1q   = 0
    j2q   = 0
    bond  = 0
    iiv   = 0
    xset  = zeros(Int8, MAX_V)

    # ── Inner helpers (closures share local state) ──────────────────────
    "Build xset[1..n]; return false if not canonical (caller → label 32)."
    function _build_xset!()
        xset[1] = Int8(1)
        jj_x = Int8(1)
        @inbounds for i in 2:n
            if state.vdeg[i - 1] != state.vdeg[i]
                jj_x += Int8(1)
            elseif state.xn[i - 1] != state.xn[i]
                jj_x += Int8(1)
            else
                ii_x = Int(state.xg[i - 1, i - 1]) - Int(state.xg[i, i])
                if ii_x > 0
                    jj_x += Int8(1)
                elseif ii_x < 0
                    return false   # not canonical
                end
            end
            xset[i] = jj_x
        end
        return true
    end

    @label dsum_outer
        # qg21:12765-12781 — advance dsum, distribute xg(i,i) greedily.
        dsum += 1
        if dsum > nloop
            return nothing
        end
        iiv = 2 * dsum
        @inbounds for i in n:-1:rhop1
            jj = min(iiv, Int(dta_arr[i]))
            state.xg[i, i] = Int8(jj)
            iiv -= jj
        end
        if iiv != 0
            # Can't fit dsum self-loops with these dta bounds — try next dsum.
            @goto dsum_outer
        end

    @label ds_init
        # qg21:12815-12818 — ds[limin, i] = vdeg[i] - xn[i] - xg[i,i].
        @inbounds for i in limin:n
            state.ds[limin, i] = state.vdeg[i] - state.xn[i] - state.xg[i, i]
        end

        # qg21:12819-12840 — build xset (canonical class index).
        if !_build_xset!()
            @goto xg_diag_backtrack
        end

        # qg21:12841 — lps[limin] = loop - dsum + 1.
        state.lps[limin] = Int8(nloop - dsum + 1)
        lin = limin

    @label row_fill
        # qg21:12842-12853 — greedy row fill at row lin.
        iiv  = Int(state.ds[lin, lin])
        bond = min(Int(str_arr[lin]), Int(state.lps[lin]))
        @inbounds for c in n:-1:(lin + 1)
            jj = min(iiv, bond, Int(state.ds[lin, c]))
            state.xg[lin, c] = Int8(jj)
            iiv -= jj
        end
        if iiv > 0
            @goto row_decrement   # label 15
        end
        @goto row_check           # label 28

    @label row_check
        # qg21:12897-12958 — row-completion check + cross-row/col canonical.
        if lin == n
            @goto emit
        end
        msum = 0
        @inbounds for i in (lin + 1):n
            ii_m = Int(state.xg[lin, i]) - 1
            if ii_m > 0
                msum += ii_m
            end
        end
        if msum >= Int(state.lps[lin])
            @goto row_decrement
        end

        # qg21:12911-12933 — cross-row canonical (only when xset[lin]==xset[lin-1]).
        if lin > limin
            if xset[lin] == xset[lin - 1]
                same_class = true
                # Check rows: i in limin..lin-2 — if xg[i,lin-1]>xg[i,lin], canonical.
                @inbounds for i in limin:(lin - 2)
                    iiv = Int(state.xg[i, lin - 1]) - Int(state.xg[i, lin])
                    if iiv > 0
                        same_class = false   # canonical → skip the strict check
                        break
                    elseif iiv < 0
                        # Should not occur in well-formed enumeration.
                        error("qg21_7 — invariant violation at row $lin")
                    end
                end
                if same_class
                    @inbounds for c in (lin + 1):n
                        iiv = Int(state.xg[lin - 1, c]) - Int(state.xg[lin, c])
                        if iiv < 0
                            col = c
                            @goto col_aux_38   # NOT canonical
                        elseif iiv > 0
                            break              # canonical (skip rest)
                        end
                    end
                end
            end
        end

        # qg21:12935-12946 — cross-col canonical for adjacent equivalent cols.
        @inbounds for c in (lin + 2):n
            if xset[c] == xset[c - 1]
                broke = false
                for i in limin:lin
                    iiv = Int(state.xg[i, c - 1]) - Int(state.xg[i, c])
                    if iiv < 0
                        col = c
                        @goto col_aux_38
                    elseif iiv > 0
                        broke = true
                        break
                    end
                end
            end
        end

        # qg21:12947-12958 — propagate ds and advance to next row.
        @inbounds for i in (lin + 1):n
            new_ds = Int(state.ds[lin, i]) - Int(state.xg[lin, i])
            new_ds < 0 && error("qg21_8 — negative ds at lin=$lin, i=$i")
            state.ds[lin + 1, i] = Int8(new_ds)
        end
        lin += 1
        state.lps[lin] = Int8(Int(state.lps[lin - 1]) - msum)
        @goto row_fill

    @label emit
        callback(state)
        @goto row_decrement   # try next topology

    # qg21:12891-12896 — entry point from cross-row/col canonical violation.
    @label col_aux_38
        aux = -1
        @inbounds for i in col:n
            aux += Int(state.xg[lin, i])
        end
        @goto col_increment   # label 23

    @label row_decrement
        # qg21:12857-12861 — backtrack one row.
        if lin == limin
            @goto xg_diag_backtrack
        end
        lin -= 1

    # qg21:12862-12890 — try to increment some column entry within row lin.
    @label row_recheck
        @inbounds for c in n:-1:(lin + 1)
            aux = Int(state.xg[lin, c]) - 1
            if aux >= 0
                col = c
                @goto col_increment
            end
        end
        @goto row_decrement

    @label col_increment
        # qg21:12877-12890 — find a column to increment, redistribute right.
        bond = min(Int(str_arr[lin]), Int(state.lps[lin]))
        found_inc = false
        @inbounds for i in (col - 1):-1:(lin + 1)
            if min(Int(state.ds[lin, i]), bond) > Int(state.xg[lin, i])
                state.xg[lin, i] += Int8(1)
                # Refill cols i+1..n with `aux` remaining.
                for k in n:-1:(i + 1)
                    jj = min(aux, bond, Int(state.ds[lin, k]))
                    state.xg[lin, k] = Int8(jj)
                    aux -= jj
                end
                found_inc = true
                break
            else
                aux += Int(state.xg[lin, i])
            end
        end
        if found_inc
            @goto row_check
        else
            @goto row_decrement
        end

    @label xg_diag_backtrack
        # qg21:12782-12798 — find rightmost xg(i,i) > 0 to backtrack.
        j1q = 0
        @inbounds for i in n:-1:rhop1
            if state.xg[i, i] > Int8(0)
                j1q = i
                break
            end
        end
        if j1q == 0
            @goto dsum_outer
        end

        # qg21:12790-12798 — find smaller-index slot to bump by 2.
        j2q = 0
        @inbounds for i in (j1q - 1):-1:rhop1
            if state.xg[i, i] < dta_arr[i]
                state.xg[i, i] += Int8(2)
                j2q = i
                break
            end
        end
        if j2q == 0
            @goto dsum_outer
        end

        # qg21:12799-12814 — redistribute lower diagonal slots.
        iiv = -2
        @inbounds for i in (j2q + 1):j1q
            iiv += Int(state.xg[i, i])
        end
        @inbounds for i in n:-1:(j2q + 1)
            jj = min(iiv, Int(dta_arr[i]))
            state.xg[i, i] = Int8(jj)
            iiv -= jj
        end
        iiv == 0 || error("qg21_3 — invalid xg-diagonal redistribution")
        @goto ds_init
end
