# Phase 18b-1: Multi-orbit Burnside summation for tree-level |M|².
#
# Called with `bundles` = one canonical representative per orbit (filtered
# by `is_emission_canonical` in solve_tree_pipeline) and `weights` = all 1.
# Under that calling convention the double-sum Σ_{i,j} w_i · w_j · T_ij
# gives each (orbit_i, orbit_j) pair exactly once, matching the handbuilt
# pattern |M|² = Σ_channels |M_c|² + cross terms (test_bhabha.jl:
# T_tt + T_ss − 2·T_int for identical-fermion exchange).
#
# The function is also correct under the full-Burnside calling convention
# (all orbit-members with w_i = |Stab_i|/|G| = 1/|Orbit_i|): the double-sum
# collapses by the orbit-stabiliser theorem to the same per-orbit-pair
# contributions.  In practice we avoid that path because cross-bundle
# spin_sum_interference (interference.jl:102) keys on bar-momentum names
# and fails when orbit-members have automorphic momentum relabelings.
#
# Diagonal term : spin_sum_amplitude_squared(line_L, line_R) — product of
#   two per-line traces (src/v2/spin_sum.jl:64-73).
# Cross term (i ≠ j) : spin_sum_interference((L_i,R_i), (L_j,R_j)) — single
#   closed-loop trace (src/v2/interference.jl:43-84).
# Fermion relative sign flows from AmplitudeBundle.fermion_sign
# (qdis_fermion_sign, qgen.jl Phase 13) multiplied per pair.  Identical-
# fermion Bhabha gets sign_s · sign_t = −1 on the s×t cross term,
# reproducing the handbuilt M = M_t − M_s anticommutation.
#
# Phase 18b-1 scope (Option A): trace-only AlgSum.  Propagator
# denominators live in `bundle.denoms`; the caller applies 1/denom at
# evaluation time.  Symbolic 1/pair(q,q) support is tracked in bead
# feynfeld-rj1l (Option B, best-in-class retrofit).

"""
    combine_m_squared_burnside(bundles, weights) -> AlgSum

Burnside-weighted |M|² (trace only) summed over all emission pairs.

Arguments:
  - `bundles` : Vector of `AmplitudeBundle` — one per emission from
    `_foreach_emission` (may include multiple members of the same orbit).
  - `weights` : Vector of `Rational{Int}` — `weight[i] = |Stab_i|/|G|`.

Returns an `AlgSum` equal to Σ_{i,j} s_i·s_j·w_i·w_j · T_ij, where
`s_i = bundles[i].fermion_sign` and T_ij is the diagonal (i = j) or
interference (i ≠ j) trace.  Denominators are NOT included — see file
header and bead feynfeld-rj1l.
"""
function combine_m_squared_burnside(bundles::Vector{AmplitudeBundle},
                                    weights::Vector{Rational{Int}})
    length(bundles) == length(weights) ||
        error("combine_m_squared_burnside: length(bundles)=$(length(bundles)) " *
              "≠ length(weights)=$(length(weights))")
    m_sq = Main.FeynfeldX.AlgSum()
    for i in 1:length(bundles), j in 1:length(bundles)
        bi, bj = bundles[i], bundles[j]
        T_ij   = _pair_trace(bi, bj, i == j)
        w      = weights[i] * weights[j]
        s      = Rational{Int}(bi.fermion_sign * bj.fermion_sign)
        m_sq   = m_sq + s * w * T_ij
    end
    m_sq
end

# Diagonal (same bundle): two-line product of per-line traces.
# Off-diagonal (different bundles): single closed-loop interference trace.
# Empty line_chains (φ³ scalar): amplitude is alg(1); every pair is alg(1).
# Phase 18b-1 supports 0- and 2-line bundles only; multi-vertex lines
# (3+ chains) are deferred to Phase 18b-3.
function _pair_trace(bi::AmplitudeBundle, bj::AmplitudeBundle, is_diagonal::Bool)
    n_i, n_j = length(bi.line_chains), length(bj.line_chains)
    if n_i == 0 && n_j == 0
        return Main.FeynfeldX.alg(1)
    end
    (n_i == 2 && n_j == 2) ||
        error("_pair_trace: Phase 18b-1 supports 0- and 2-line bundles only " *
              "(got $n_i × $n_j lines); multi-vertex fermion lines deferred " *
              "to Phase 18b-3")
    if is_diagonal
        Main.FeynfeldX.spin_sum_amplitude_squared(bi.line_chains[1], bi.line_chains[2])
    else
        Main.FeynfeldX.spin_sum_interference((bi.line_chains[1], bi.line_chains[2]),
                                             (bj.line_chains[1], bj.line_chains[2]))
    end
end
