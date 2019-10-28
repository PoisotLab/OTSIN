using Distributions: MultivariateNormal


function optimaltransport_safe(M::AbstractMatrix, a::AbstractVector, b::AbstractVector;
                        λ::Real=1.0, maxitter::Int=10, ϵ::Real=1e-5)
    Q = exp.(λ * M)
    iter = 0
    while iter < maxitter && maximum(abs.(sum(Q, dims=2)[:] - a)) > ϵ
        iter += 1
        Q = a .* (sum(Q, dims=2) .\ Q)
        Q = Q .* transpose(b) ./ sum(Q, dims=1)
    end
    return Q
end

function laplace_posterior(M, P; γ=1e-2,
                maxitter::Int=20, tol=1e-4)
    a, b = marginals(P)
    Q = M -> Q = optimaltransport_safe(exp.(M), a, b)
    l = M -> relative_entropy(P, Q(M)) + γ * sum(M .* M)
    hess = Zygote.hessian(l, M)
    return MultivariateNormal(vec(exp.(M)), Symmetric(hess+tol*I)^-1)
end

#TODO: test and explore this.
