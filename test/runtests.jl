using Test
using Feynfeld

@testset "Feynfeld.jl" begin
    @testset "Algebra" begin
        include("algebra/test_types.jl")
    end
end
