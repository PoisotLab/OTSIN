#=
Created on August 2019
Last update: Tuesday 01 October 2019

@author: Michiel Stock
michielfmstock@gmail.com

Methods for fitting a utility matrix to an observed interaction network.
=#


using Zygote, Optim


# regularization methods
abstract type Regularization end

struct L1 <: Regularization end
struct L2 <: Regularization end
struct ElasticNet <: Regularization
    α::Real
    function ElasticNet(α::Real=0.5)
        @assert α ≥ 0 && α ≤ 1
        return new(α)
    end
end
struct Exp <: Regularization end
struct Entropic <: Regularization
    offset::Real
end

Entropic() = Entropic(0.0)

r(x, reg::L1) = sum(abs.(x))
r(x, reg::L2) = sum(x.^2)
r(x, reg::ElasticNet) = reg.α * sum(abs.(x)) + (1 - reg.α) * sum(x.^2)
r(x, reg::Exp) = sum(exp.(x))
r(x, reg::Entropic) = sum(x .* log.(x .+ reg.offset))


function _fitM(P::AbstractMatrix; γ=1e-2, reg::Regularization=L2())

    l = M -> relative_entropy(P, optimaltransport(M)) + γ * r(M, reg)
    dM = similar(P)
    ∇l!(dM, M) = dM .= l'(M)
    return Optim.minimizer(optimize(l, ∇l!, zero(P), BFGS()))
end

function _fitM(P::AbstractMatrix, a::AbstractVector; γ=1e-2, reg::Regularization=L2())

    l = M -> relative_entropy(P, optimaltransport(M, a)) + γ * r(M, reg)
    dM = similar(P)
    ∇l!(dM, M) = dM .= l'(M)
    return Optim.minimizer(optimize(l, ∇l!, zero(P), BFGS()))
end

function _fitM(P, a::AbstractVector, b::AbstractVector; γ=1e-2,
                                    reg::Regularization=L2(),
                                    maxitter::Int=100, ϵ::Real=1e-5)

    l = M -> relative_entropy(P, optimaltransport(M, a, b; maxitter=maxitter,
                                        ϵ=ϵ)) + γ * r(M, reg)
    dM = similar(P)
    ∇l!(dM, M) = dM .= l'(M)
    return Optim.minimizer(optimize(l, ∇l!, zero(P), BFGS()))
end

function _fitMnograd(P, a::AbstractVector, b::AbstractVector; γ=1e-2,
                                    reg::Regularization=L2(),
                                    maxitter::Int=100, ϵ::Real=1e-5)

    l = M -> relative_entropy(P, optimaltransport(M, a, b; maxitter=maxitter,
                                        ϵ=ϵ)) + γ * r(M, reg)
    dM = similar(P)
    return Optim.minimizer(optimize(l, zero(P), BFGS()))
end

"""
    fitM(P::AbstractMatrix; fix_a::Bool, fix_b::Bool,
                    γ=1e-2,
                    reg::Regularization=L2(),
                    maxitter::Int=100, ϵ::Real=1e-5)

Fit a utility matrix `M` to an observed interaction frequency matrix `P`.
The parameters `fix_a` and `fix_b` determine whether the marginals should be
fixed.
"""
function fitM(P::AbstractMatrix; fix_a::Bool=true, fix_b::Bool=true,
                    γ=1e-2,
                    reg::Regularization=L2(),
                    maxitter::Int=100, ϵ::Real=1e-5)
    a, b = marginals(P)
    if fix_a && fix_b
        return _fitM(P, a, b, reg=reg, γ=γ, maxitter=maxitter, ϵ=ϵ)
    elseif fix_a
        return _fitM(P, a, reg=reg, γ=γ)
    elseif fix_b
        return _fitM(P', b, reg=reg, γ=γ)'
    else
        return _fitM(P, reg=reg, γ=γ)
    end
end

# loss over several observed frequency matrices
function loss(M, Ps::Array{M, 1} where M <: AbstractMatrix,
                fix_a::Bool, fix_b::Bool;
                γ=1e-2,
                reg::Regularization=L2(),
                maxitter::Int=10,
                ϵ::Real=1e-5)
    l = 0.0
    for P in Ps
        a, b = marginals(P)
        Q = optimaltransport(M, fix_a ? a : nothing,
                                fix_b ? b : nothing,
                                ϵ=ϵ, maxitter=maxitter)
        l += relative_entropy(P, Q)
    end
    l += γ * r(M, reg)
    return l
end


"""
    fitM(Ps::Array{M, 1} where M <: AbstractMatrix,
                    fix_a::Bool, fix_b::Bool,
                    γ=1e-2,
                    reg::Regularization=L2(),
                    maxitter::Int=100, ϵ::Real=1e-5)

Fit a utility matrix `M` to a collection of observed probability matrices `Ps=[Pi]`.
The parameters `fix_a` and `fix_b` determine whether the marginals should be
fixed.
"""
function fitM(Ps::Array{M, 1} where M <: AbstractMatrix;
                fix_a::Bool, fix_b::Bool,
                γ=1e-2,
                reg::Regularization=L2(),
                maxitter::Int=20,
                ϵ::Real=1e-4, usegrad=true)
    l = M -> loss(M, Ps, fix_a, fix_b, γ=γ, reg=reg,
                        maxitter=maxitter, ϵ=ϵ)
    dM = similar(first(Ps))
    if usegrad
        ∇l!(dM, M) = dM .= l'(M)
        return Optim.minimizer(optimize(l, ∇l!, zero(first(Ps)), BFGS()))
    else
        return Optim.minimizer(optimize(l, zero(first(Ps)), BFGS()))
    end
end

using StatsBase: sample

"""
    fitMsgd(Ps::Array{M, 1} where M <: AbstractMatrix;
                    fix_a::Bool, fix_b::Bool,
                    γ=1e-2,
                    reg::Regularization=L2(),
                    maxitter::Int=20,
                    ϵ::Real=1e-4,
                    α::Real=0.1,
                    β::Real=0.8,
                    nsteps::Real=100)

Fit a utility matrix `M` to a collection of observed probability matrices
`Ps=[Pi]` using stochastic gradient descent with momentum.
The parameters `fix_a` and `fix_b` determine whether the marginals should be
fixed.
"""
function fitMsgd(Ps::Array{M, 1} where M <: AbstractMatrix;
                fix_a::Bool, fix_b::Bool,
                γ=1e-2,
                reg::Regularization=L2(),
                maxitter::Int=20,
                ϵ::Real=1e-4,
                α::Real=0.1,
                β::Real=0.8,
                nsteps::Real=100)  # momentum parameter)
    global_loss = M -> loss(M, Ps, fix_a, fix_b, γ=γ, reg=reg,
                                    maxitter=maxitter, ϵ=ϵ)
    losses = zeros(nsteps)
    ΔM = similar(first(Ps))
    fill!(ΔM, 0.0)
    M = zeros(size(ΔM)...)
    for i in 1:nsteps
        losses[i] = global_loss(M)
        for rep in 1:length(Ps)
            P = sample(Ps)
            a, b = marginals(P)
            l = M -> relative_entropy(P,
                    optimaltransport(M, fix_a ? a : nothing,
                            fix_b ? b : nothing,
                            maxitter=maxitter, ϵ=ϵ)) + (γ / length(Ps)) * r(M, reg)
            ΔM .*= 1 - β
            ΔM .+= α * l'(M)
            M .-= ΔM
        end
    end
    return M, losses
end
