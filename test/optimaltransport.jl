@testset "optimal transport" begin

    P = randn(4, 5).^2 |> x -> x .= x / sum(x)

    M = 0.1randn(size(P))

    a, b = marginals(P)

    @test isprobability(P)
    @test !isprobability(P .+ 0.1)
    @test isprobability(a)
    @test isprobability(b)

    @test entropy(P[1,:] |> OTSIN.normalize) ≈ entropy(P[1,:], norm=true)

    @testset "free OT" begin
        Q = optimaltransport(M, λ=1.0)
        @test isprobability(Q)
        @test entropy(Q) > entropy(optimaltransport(M, λ=10.0))
        @test relative_entropy(P, Q) > 0
        @test Q ≈ optimaltransport(M, nothing, nothing, λ=1.0)
    end

    @testset "one marginal OT" begin
        Q = optimaltransport(M, a, λ=1.0)
        @test isprobability(Q)
        @test marginals(Q)[1] ≈ a
        @test Q ≈ optimaltransport(M, a, nothing, λ=1.0)
    end

    @testset "Sinkhorn" begin
        Q = optimaltransport(M, a, b, λ=1.0, ϵ=1e-10)
        @test isprobability(Q)
        @test all(marginals(P) .≈ marginals(Q))
        @test optimaltransport(M, a, b, λ=0) ≈ optimaltransport(a, b)
        @test Q ≈ optimaltransport(M, a, b, λ=1.0, ϵ=1e-10)
    end

    #  check if some elements of a and b can be zero
    @testset "Sinkhorn sparse" begin
        Q = optimaltransport(M, [0.5, 0.2, 0.0, 0.3],
                    [0.1, 0.1, 0.1, 0.0, 0.7], λ=1.0, ϵ=1e-10)
        a_, b_ = marginals(Q)
        @assert a_[3] < 1e-6
        @assert b_[4] < 1e-6
        @assert b_[1] ≈ 0.1
        @assert sum(Q) ≈ 1.0
    end

    @testset "perturbation" begin
        ∇af, ∇bf = perturbation(Q -> sum(Q.^2), M, a, b)
        @test ∇af ⋅ ones(size(∇af)) < 1e-6
        @test ∇bf ⋅ ones(size(∇bf)) < 1e-6
    end

end
