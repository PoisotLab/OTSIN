
P = rand(5, 6) |> OTSIN.normalize
Ps = [rand(size(P)...) |> OTSIN.normalize for i in 1:10]
M = randn(size(P)...)

@testset "fitting" begin

    @testset "regularization" begin
        @test r(M, L2()) ≈ r(M, ElasticNet(0.0)) < r(10M, L2())
        @test r(M, L1()) ≈ r(M, ElasticNet(1.0)) < r(10M, L1())
        @test r(M, ElasticNet(0.5)) ≈ 0.5r(M, L1()) + 0.5r(M, L2())
    end

    @testset "fitting" begin
        Mfit = fitM(P, fix_a=false, fix_b=true, γ=0.01)
        @test sum(Mfit) < 1e-5  # approximately 0
        Mfit = fitM(Ps, fix_a=true, fix_b=true)
        Mfit, losses = fitMsgd(Ps, fix_a=true, fix_b=true)
    end

    @testset "uneven marginals" begin
        # a is fixed, b not
        a = rand(10) |> OTSIN.normalize
        Ps = [a * (OTSIN.normalize(rand(1, 5))) for i in 1:10]
        Mfree, losses = fitMsgd(Ps, fix_a=false, fix_b=false)
    end
end
