
using Distributions: MultivariateNormal

function optimaltransport_(M::AbstractMatrix, a::AbstractVector, b::AbstractVector;
                        λ::Real=1.0, maxitter::Int=100, ϵ::Real=1e-5)
    Q = exp.(λ * M)
    iter = 0
    while iter < maxitter && maximum(abs.(sum(Q, dims=2)[:] - a)) > ϵ
        iter += 1
        Q = Diagonal(a ./ sum(Q, dims=2)[:]) * Q
        Q = Q * Diagonal(b ./ sum(Q, dims=1))
    end
    return Q
end


function laplace_posterior(M̂, P, γ=1e-2,
                reg::Regularization=L2(),
                maxitter::Int=20,
                ϵ::Real=1e-4)
    a, b = marginals(P)
    Q = M -> Q = optimaltransport_(M, a, b,
                    ϵ=ϵ, maxitter=maxitter)
    l = M -> KL(P, Q(M)) + γ * r(M, reg)
    hess = Zygote.hessian(l, M̂)
    return MultivariateNormal(vec(M̂), hess^-1)
end

#TODO: test and explore this.
