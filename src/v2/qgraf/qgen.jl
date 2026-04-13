#  qgen — field assignment via dpntro lookup.
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:13797-14464
#
#  Phase 12a — dpntro builder (this file).
#  Phase 12b — vmap/lmap construction (extension to qg10).
#  Phase 12c — qgen recursive backtracker.

"""
    build_dpntro(vertex_rules) -> Dict{Int, Vector{Vector{Symbol}}}

Build the qgen lookup table.  Input is any iterable of `Vector{Symbol}`
representing expanded vertex rules (each Vector already sorted by
`_expand_vertex`'s `sort` step).  Output buckets the rules by arity and
sorts each bucket lexicographically for deterministic iteration.

Source: qgraf-4.0.6.f08:13889
  "vfo(vv) = stib(stib(dpntro(0)+vdeg(vv))+pmap(vv,1))"

The qgraf two-level nesting (degree → first-particle → list) is collapsed
to a single Dict-of-lists; downstream qgen filters at lookup time.  A
two-level dict is a future optimisation when profiling shows it matters.
"""
function build_dpntro(vertex_rules)
    dp = Dict{Int, Vector{Vector{Symbol}}}()
    for rule in vertex_rules
        d = length(rule)
        push!(get!(Vector{Vector{Symbol}}, dp, d), Vector{Symbol}(rule))
    end
    for d in keys(dp)
        sort!(dp[d])
    end
    dp
end
