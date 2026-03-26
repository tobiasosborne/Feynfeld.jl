# Tests for Minkowski registry setup (TensorGR bridge)
# Validates metric contraction, dimension traces, and TensorGR integration.

@testset "Minkowski registry (4D)" begin
    reg = minkowski_registry(dim=4)
    @test reg isa TensorGR.TensorRegistry
    @test TensorGR.has_manifold(reg, :M4)
    @test TensorGR.has_tensor(reg, :η)
    @test TensorGR.has_tensor(reg, :δ)

    # Check manifold dimension
    mp = TensorGR.get_manifold(reg, :M4)
    @test mp.dim == 4

    # Metric is flat
    @test TensorGR.is_flat(reg, :η)

    # Metric contraction: η^μν η_νρ = δ^μ_ρ
    TensorGR.with_registry(reg) do
        η_up = TensorGR.Tensor(:η, [TensorGR.up(:μ), TensorGR.up(:ν)])
        η_dn = TensorGR.Tensor(:η, [TensorGR.down(:ν), TensorGR.down(:ρ)])
        expr = TensorGR.tproduct(1 // 1, TensorGR.TensorExpr[η_up, η_dn])
        result = TensorGR.contract_metrics(expr)
        # Should reduce to δ^μ_ρ
        @test result isa TensorGR.Tensor
        @test result.name == :δ
    end

    # Metric trace: η^μ_μ = 4
    TensorGR.with_registry(reg) do
        η_mixed = TensorGR.Tensor(:η, [TensorGR.up(:μ), TensorGR.down(:μ)])
        result = TensorGR.contract_metrics(η_mixed)
        @test result isa TensorGR.TScalar
        @test result.val == 4 // 1
    end
end

@testset "Minkowski registry (D-dimensional)" begin
    reg = minkowski_registry(dim=:D)
    mp = TensorGR.get_manifold(reg, :M4)
    @test mp.dim === :D

    # Metric trace: η^μ_μ = D (symbolic)
    TensorGR.with_registry(reg) do
        η_mixed = TensorGR.Tensor(:η, [TensorGR.up(:μ), TensorGR.down(:μ)])
        result = TensorGR.contract_metrics(η_mixed)
        @test result isa TensorGR.TScalar
        @test result.val === :D
    end
end
