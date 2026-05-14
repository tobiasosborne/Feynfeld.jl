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
    m_sq = AlgSum()
    for i in 1:length(bundles), j in 1:length(bundles)
        bi, bj = bundles[i], bundles[j]
        T_ij   = _pair_trace(bi, bj, i == j)
        w      = weights[i] * weights[j]
        s      = Rational{Int}(bi.fermion_sign * bj.fermion_sign)
        add!(m_sq, T_ij, s * w)
    end

    # Phase 18b-4: external-boson polarisation sums (Feynman gauge,
    # Σ_λ ε^μ ε^{ν*} = -g^{μν}). Every bundle of a given process carries
    # the same canonical `:eps_<leg>` indices (see emission_to_amplitude),
    # so the sum factorises out of the i,j double sum and is applied once.
    # The conjugate amplitude relabels each index `x → x_` (spin_sum.jl
    # `_conjugate_gammas`), so the metric ties `:eps_<leg>` to
    # `:eps_<leg>_`; `solve_tree_pipeline`'s `contract` does the rest.
    if !isempty(bundles)
        pol_idxs = bundles[1].boson_pols
        for b in bundles
            b.boson_pols == pol_idxs ||
                error("combine_m_squared_burnside: bundles disagree on external " *
                      "boson polarisation indices ($(b.boson_pols) vs $pol_idxs)")
        end
        for mu in pol_idxs
            mu_conj = LorentzIndex(Symbol(mu.name, :_), mu.dim)
            m_sq = m_sq * polarization_sum(mu, mu_conj)
        end
    end
    m_sq
end

# Per-(i,j) trace, dispatched on fermion-line count:
#   0×0  (φ³ scalar)      : fermion part is alg(1).
#   1×1  (Compton, qq̄→gg quark line) : one external-spinor pair, boson
#         polarisation indices free. Diagonal is Σ_spins|M_i|²; the
#         off-diagonal is the 1-line interference Σ_spins M_j* M_i —
#         valid because every emission of a 1-line process pins the same
#         external bar/plain momenta (Phase 18b-4).
#   2×2  (ee→μμ, Bhabha)  : product of two per-line traces (diagonal) or
#         the reconnected closed-loop interference trace (off-diagonal).
# Mixed and 3+-line bundles are not yet supported.
#
# The boson sub-amplitude (off-fermion-line vertex factors, e.g. the
# triple-gluon vertex in qq̄→gg s-channel) multiplies the fermion trace:
# M_j contributes `boson_factor_j`, M_i* the conjugate (every Lorentz
# index relabelled x → x_, matching the fermion line's `_conjugate_gammas`).
# For fermion-only bundles both factors are alg(1) and this is inert.
function _pair_trace(bi::AmplitudeBundle, bj::AmplitudeBundle, is_diagonal::Bool)
    n_i, n_j = length(bi.line_chains), length(bj.line_chains)
    fermion = if n_i == 0 && n_j == 0
        alg(1)
    elseif n_i == 1 && n_j == 1
        is_diagonal ?
            _single_line_trace(bi.line_chains[1]) :
            _line_trace(bj.line_chains[1], bi.line_chains[1])
    elseif n_i == 2 && n_j == 2
        is_diagonal ?
            spin_sum_amplitude_squared(bi.line_chains[1], bi.line_chains[2]) :
            spin_sum_interference((bi.line_chains[1], bi.line_chains[2]),
                                  (bj.line_chains[1], bj.line_chains[2]))
    else
        error("_pair_trace: unsupported fermion-line counts $n_i × $n_j " *
              "(supported: 0×0, 1×1, 2×2); mixed / 3+-line bundles not yet handled")
    end
    fermion * bj.boson_factor * _conjugate_algsum_indices(bi.boson_factor)
end
