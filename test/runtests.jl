using Test
using Feynfeld
import TensorGR

@testset "Feynfeld.jl" begin
    @testset "Dimensions" begin
        include("algebra/test_dimensions.jl")
    end
    @testset "Algebra" begin
        include("algebra/test_types.jl")
    end
    @testset "Momentum" begin
        include("algebra/test_momentum.jl")
    end
    @testset "Pair" begin
        include("algebra/test_pair.jl")
    end
    @testset "SPContext" begin
        include("algebra/test_sp_context.jl")
    end
    @testset "ExpandScalarProduct" begin
        include("algebra/test_expand_sp.jl")
    end
    @testset "Eps" begin
        include("algebra/test_eps.jl")
    end
    @testset "Contract" begin
        include("algebra/test_contract.jl")
    end
    @testset "MUnit Lorentz" begin
        include("algebra/test_munit_lorentz.jl")
    end
    @testset "Minkowski" begin
        include("algebra/test_minkowski.jl")
    end
end
