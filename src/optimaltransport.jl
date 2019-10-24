#=
Created on 13 August 2019
Last update: 30 September 2019

@author: Michiel Stock
michielfmstock@gmail.com

Functions to simulate species interaction networks with optimal transport.
=#


"""
    zerodiv(a, b)

Division which results 0 if both arguments are 0.
"""
zerodiv(a, b) = a == b == zero(a) ? zero(a) : a / b

function optimaltransport_(M::AbstractMatrix, a::AbstractVector; λ::Real=1.0)
    @assert length(a) == size(M, 1) throw(DimensionMismatch(
                "a ($(size(a))) does not match M ($(size(M)))"))
    K = exp.(λ * M)
    return Diagonal(a ./ sum(K, dims=2)[:]) * K
end


function optimaltransport_inplace(M::AbstractMatrix, a::AbstractVector, b::AbstractVector;
                        λ::Real=1.0, maxitter::Int=100, ϵ::Real=1e-5)
    Q = exp.(λ * M)
    iter = 0
    while iter < maxitter && maximum(abs.(sum(Q, dims=1)[:] - b)) > ϵ
        iter += 1
        lmul!(Diagonal(a ./ sum(Q, dims=2)[:]), Q)
        rmul!(Q, Diagonal(b ./ sum(Q, dims=1)[:]))
    end
    return Q
end

# version that is differentiable
function optimaltransport_(M::AbstractMatrix, a::AbstractVector, b::AbstractVector;
                        λ::Real=1.0, maxitter::Int=100, ϵ::Real=1e-5)
    Q = exp.(λ * M)
    iter = 0
    while iter < maxitter && maximum(abs.(sum(Q, dims=2)[:] - a)) > ϵ
        iter += 1
        Q = Diagonal(zerodiv.(a, sum(Q, dims=2)[:])) * Q
        Q = Q * Diagonal(zerodiv.(b, sum(Q, dims=1)))
    end
    return Q
end

function optimaltransport(M::AbstractMatrix,
                a::Union{AbstractVector,Nothing}=nothing,
                b::Union{AbstractVector,Nothing}=nothing;
                λ::Real=1.0, maxitter::Int=1000, ϵ::Real=1e-5)
    if a == nothing && b == nothing
        return exp.(λ * M) |> x -> x / sum(x)
    elseif b == nothing
        return optimaltransport_(M, a, λ=λ)
    elseif a == nothing
        return optimaltransport_(M', b, λ=λ)'
    else
        return optimaltransport_(M, a, b, λ=λ, maxitter=maxitter, ϵ=ϵ)
    end
end

optimaltransport(a::AbstractVector, b::AbstractVector; λ::Real=0.0) = a * b'
