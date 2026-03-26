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
    @testset "Minkowski" begin
        include("algebra/test_minkowski.jl")
    end
end
