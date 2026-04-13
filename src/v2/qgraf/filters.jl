#  qgraf filter predicates.
#
#  Source map (qgraf-4.0.6.f08):
#    qumpi(1,...)  3690  → onepi / nobridge
#    qumpi(2,...)         → nosbridge
#    qumpi(3,...)         → notadpole
#    qumpi(4,...)         → onshell
#    qumvi(1,...)  3777  → nosnail
#    qumvi(2,...)         → onevi
#    qumvi(3,...)         → onshellx
#    qgsig         13669 → nosigma
#    qcyc          18830 → cycli
#    inline in qg21       → noselfloop, nodiloop, noparallel
#
#  Phase 14a: the three inline xg-pattern filters (smallest, all in this file).

"""
    has_no_selfloop(state) -> Bool

True iff no internal vertex carries a self-loop (xg[i,i] == 0 for all
internal i).  Source: qg21:13960-13967 (dflag(dflaty%nosl) check).
"""
function has_no_selfloop(state::TopoState)
    @inbounds for i in Int(state.rhop1):Int(state.n)
        state.xg[i, i] > Int8(0) && return false
    end
    return true
end

"""
    has_no_diloop(state) -> Bool

True iff no two distinct vertices share more than one edge between them
(xg[i,j] ≤ 1 for all i < j).  Source: qg21:13969-13978 (dflag(dflaty%nodl)).
"""
function has_no_diloop(state::TopoState)
    n = Int(state.n)
    @inbounds for i in 1:(n - 1)
        for j in (i + 1):n
            state.xg[i, j] > Int8(1) && return false
        end
    end
    return true
end

"""
    has_no_parallel(state) -> Bool

True iff `has_no_diloop` AND no internal vertex has xg[i,i] > 3 (more
than one self-loop pair).  Source: qg21:13065-13076 (dflag(dflaty%nopa)).
"""
function has_no_parallel(state::TopoState)
    has_no_diloop(state) || return false
    @inbounds for i in Int(state.rhop1):Int(state.n)
        state.xg[i, i] > Int8(3) && return false
    end
    return true
end
