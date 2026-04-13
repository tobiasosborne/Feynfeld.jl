#  Phase 15: symmetry factor — S_local component.
#  Source: qgraf-4.0.6.f08:14361-14411 (sets ndsym(symt%l)).
#  Cross-ref: Nogueira J.Comp.Phys 105 (1993) 279, p. 281 §3.

using Test
include("../../../src/v2/FeynfeldX.jl")
using .FeynfeldX.QgrafPort: Partition, TopoState, MAX_V,
                              qg21_enumerate!, compute_qg10_labels,
                              compute_local_sym_factor

@testset "Phase 15: S_local symmetry factor" begin

    # Build a pmap holding all-:phi for a phi3 topology.
    function _all_phi_pmap(state)
        n = Int(state.n)
        pmap = fill(:_, n, MAX_V)
        @inbounds for v in 1:n
            for slot in 1:Int(state.vdeg[v])
                pmap[v, slot] = :phi
            end
        end
        pmap
    end

    @testset "phi3 tree φ→φφ — S_local = 1 (no parallel edges, no self-loops)" begin
        # Single internal vertex (4) connected to all 3 externals.
        # rdeg(4) = 3, vdeg(4) = 3, j starts at 4 > vdeg → empty loop.
        p = Partition(Int8(3), Int8[1], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        @test length(emits) == 1
        labels = compute_qg10_labels(emits[1])
        pmap   = _all_phi_pmap(emits[1])
        conj   = Dict{Symbol, Symbol}(:phi => :phi, :_ => :_)
        @test compute_local_sym_factor(emits[1], labels, pmap, conj) == 1
    end

    @testset "phi3 tree φφ→φφ — S_local = 1 (single internal edge, no multi)" begin
        # Vertex 5 has rdeg=2, vdeg=3.  j=3: gam(5,6)=1 → k=4; one slot, kk=1,
        # no multiplication.  Vertex 6 similar (already-visited).  S_local=1.
        p = Partition(Int8(4), Int8[2], Int8(3), Int8(0))
        s = TopoState(p)
        emits = TopoState[]
        qg21_enumerate!(s) do state
            push!(emits, deepcopy(state))
        end
        @test length(emits) == 1
        labels = compute_qg10_labels(emits[1])
        pmap   = _all_phi_pmap(emits[1])
        conj   = Dict{Symbol, Symbol}(:phi => :phi, :_ => :_)
        @test compute_local_sym_factor(emits[1], labels, pmap, conj) == 1
    end

    # NB: vacuum (n_ext=0) self-loop cases require manual labels — qgraf's
    # qg10 (and our compute_qg10_labels) error on vaux=0 at the first pick.
    # Non-vacuum self-loop cases (e.g. 1-ext 1-loop self-energy) need
    # field-rules that allow the topology, which is a Phase 16+ concern.
    # The self-loop branch of compute_local_sym_factor is therefore exercised
    # only in integration tests once full self-loop topologies become
    # producible by the pipeline.

end
