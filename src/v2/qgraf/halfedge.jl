#  Phase 18a-2: half-edge labelling (compute_amap).
#
#  Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12133-12158
#         (internal-vertex slot loop, triple-case logic)
#         f08:12342-12344 (external back-write)
#         f08:13320-13340 (ege table — lex (i,j) order from rhop1)
#
#  amap[v, slot] = global half-edge id.  Externals 1..nleg; internals
#  nleg+1..nli in lex order of internal vertex pairs.  Both half-edges
#  of any single edge carry the same id (pairing invariant).

"""
    compute_amap(state, labels) -> Matrix{Int}

Build the half-edge labelling matrix `amap[v, slot]` per qgraf's qg10
construction.

External half-edge of vertex i (1 ≤ i ≤ nleg) has slot 1 with
`amap[i, 1] = i`; the same id is back-written to its internal-vertex
peer slot via `vmap`/`lmap` (f08:12342-12344).

Internal half-edges (vmap[v, s] > nleg) are labelled by the global edge
id, computed as `ege[min(v,j3), max(v,j3)] + offset`. The offset is:
- 0 for a single edge (gam == 1)
- (slot − 1 − rdeg[v]) ÷ 2 for a self-loop pair (qgraf's integer
  division pairs the two half-edge slots to the same id)
- the backward-scan count of consecutive `vmap[v, *] == j3` slots above,
  for parallel edges (gam > 1, v ≠ j3)

Pairing invariant: amap[v, s] == amap[vmap[v, s], lmap[v, s]] for every
valid (v, s).

Source: refs/qgraf/v4.0.6/qgraf-4.0.6.dir/qgraf-4.0.6.f08:12133-12158
"""
function compute_amap(state::TopoState, labels)
    n     = Int(state.n)
    n_ext = Int(state.n_ext)
    rhop1 = Int(state.rhop1)

    # ege[i, j] = base global edge id for vertex pair (i, j), i ≤ j,
    # both internal. Lex (i, j) order starting at rhop1 = nleg + 1.
    ege = zeros(Int, n, n)
    next_id = rhop1
    for i in rhop1:n
        for j in i:n
            mult_raw = Int(state.xg[i, j])
            mult_raw == 0 && continue
            ege[i, j] = next_id
            mult = i == j ? mult_raw ÷ 2 : mult_raw
            next_id += mult
        end
    end

    amap = zeros(Int, n, MAX_V)

    # External half-edges + back-write to peer (f08:12342-12344).
    for i in 1:n_ext
        amap[i, 1] = i
        peer_v = Int(labels.vmap[i, 1])
        peer_s = Int(labels.lmap[i, 1])
        amap[peer_v, peer_s] = i
    end

    # Internal half-edges (f08:12133-12158).
    for v in rhop1:n
        for s in 1:Int(state.vdeg[v])
            j3 = Int(labels.vmap[v, s])
            j3 > n_ext || continue                # external slot, already set
            j1, j2 = min(v, j3), max(v, j3)
            base = ege[j1, j2]
            mult = Int(state.xg[j1, j2])
            offset = if mult == 1
                0
            elseif v == j3
                (s - 1 - Int(labels.rdeg[v])) ÷ 2
            else
                kk = 0
                for s_back in (s - 1):-1:1
                    Int(labels.vmap[v, s_back]) == j3 || break
                    kk += 1
                end
                kk
            end
            amap[v, s] = base + offset
        end
    end

    amap
end
