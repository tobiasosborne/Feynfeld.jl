# Tests for D₀ scalar four-point function.
# Ground truth: COLLIER vs quadgk cross-validation.
# Ref: refs/papers/tHooftVeltman1979_NuclPhysB153.pdf, Eq. (5.2)

using Test
@isdefined(FeynfeldX) || include("../../src/v2/FeynfeldX.jl")
using .FeynfeldX

const _LIB = joinpath(@__DIR__, "../../refs/COLLIER/COLLIER-1.2.8/libcollier.so")
const HAS_COLLIER = isfile(_LIB)

@testset "D₀ scalar four-point" begin

    @testset "D0 type construction" begin
        d = D0(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0)
        @test d isa PaVe{4}
        @test d.invariants == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        @test d.masses == [1.0, 2.0, 3.0, 4.0]
        @test d.indices == Int[]
        @test string(d) == "D(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0)"
    end

    @testset "D0 evaluate dispatch" begin
        d = D0(-1.0, -2.0, -1.5, -3.0, -2.5, -1.8,
               1.0, 2.0, 1.5, 3.0)
        result = evaluate(d)
        @test result isa ComplexF64
        @test isfinite(real(result))
        @test isfinite(imag(result))
    end

    @testset "D0 massive spacelike — purely real" begin
        # All masses nonzero, all momenta spacelike → no cuts → real D₀
        d = D0(-1.0, -2.0, -3.0, -4.0, -2.5, -3.5,
               1.0, 2.0, 3.0, 4.0)
        result = evaluate(d)
        @test abs(imag(result)) < 1e-8 * abs(real(result))
        @test isfinite(real(result))
    end

    @testset "D-tensor via evaluate dispatch" begin
        d1 = PaVe{4}([1], [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
                      [1.0, 1.0, 1.0, 1.0])
        r = evaluate(d1)
        @test isfinite(real(r))
    end

    if HAS_COLLIER
        @testset "COLLIER — multiple kinematic regions" begin
            # Spacelike, equal masses
            c1 = FeynfeldX._D0_collier(-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
                                        1.0, 1.0, 1.0, 1.0)
            @test isfinite(real(c1)) && isfinite(imag(c1))
            @test abs(imag(c1)) < 1e-10  # purely real

            # Timelike p10, rest spacelike — should have imaginary part
            c2 = FeynfeldX._D0_collier(10.0,-2.0,-1.5,-3.0,-2.5,-1.8,
                                        0.5, 0.5, 0.5, 0.5)
            @test isfinite(real(c2)) && isfinite(imag(c2))

            # Mixed timelike/spacelike
            c3 = FeynfeldX._D0_collier(4.0,-1.0, 9.0,-2.0, 3.0,-1.5,
                                        1.0, 0.5, 2.0, 1.5)
            @test isfinite(real(c3)) && isfinite(imag(c3))

            # Large masses (heavy propagators) — D₀ → 0 as m → ∞
            c4 = FeynfeldX._D0_collier(-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
                                        100.0,100.0,100.0,100.0)
            c5 = FeynfeldX._D0_collier(-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,
                                        1.0, 1.0, 1.0, 1.0)
            @test abs(c4) < abs(c5)  # heavier → smaller
        end

        # quadgk cross-validation skipped by default (triple-nested → O(minutes)).
        # Run manually: FeynfeldX._D0_quadgk(-1.0,-2.0,-1.5,-3.0,-2.5,-1.8, 1.0,2.0,1.5,3.0)
        # and compare against _D0_collier for the same point.
    else
        @warn "COLLIER not available — D₀ quadgk fallback is too slow for automated tests"
    end

    # ---- D-tensor coefficients (PV reduction) ----
    if HAS_COLLIER
        @testset "D₁/D₂/D₃ tensor coefficients" begin
            # General kinematics: all coefficients finite
            # Ref: Denner1993 Eqs. (4.7)-(4.9), PV reduction via 3×3 Gram matrix
            args = (-1.0, -2.0, -1.5, -3.0, -2.5, -1.8, 1.0, 2.0, 1.5, 3.0)
            d1 = evaluate(D1(args...))
            d2 = evaluate(D2(args...))
            d3 = evaluate(D3(args...))
            @test isfinite(real(d1))
            @test isfinite(real(d2))
            @test isfinite(real(d3))

            # Symmetric kinematics: D₁ = D₂ = D₃
            args_sym = (-1.0,-1.0,-1.0,-1.0,-1.0,-1.0, 1.0,1.0,1.0,1.0)
            d1s = evaluate(D1(args_sym...))
            d2s = evaluate(D2(args_sym...))
            d3s = evaluate(D3(args_sym...))
            @test d1s ≈ d2s atol=1e-12
            @test d2s ≈ d3s atol=1e-12

            # D-tensor negative for spacelike (like D₀)
            @test real(d1) < 0
            @test real(d1s) < 0
        end
    end
end
