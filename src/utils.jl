"""
    isprobability(P)

Test if `P` is a is a probability distribution.
"""
isprobability(P) = all(P .≥ 0.0) && sum(P) ≈ 1.0

"""
    marginals(P::AbstractMatrix)

Returns the marginals of a probabilistic matrix P.
"""
marginals(P::AbstractMatrix) = sum(P, dims=2)[:], sum(P, dims=1)[:]

"""
    marginals(Ps::AbstractMatrix)

Returns the marginals of as array of probabilistic matrices Ps.
"""
function marginals(Ps::Array{M, 1} where M <: AbstractMatrix)
    A = hcat([sum(P, dims=2) for P in Ps]...)
    B = hcat([sum(P, dims=1)' for P in Ps]...)
    return A, B
end

"""
    uniform(n::Int)

Returns the PMF of length n, i.e. a vector of length n filled with 1/n.
"""
uniform(n::Int) = ones(n) / n

normalize(x) = x / sum(x)
#=
struct UtilityMatrix <: AbstractMatrix
    M::AbstractMatrix
    bottom_species::AbstractVector
    top_species::AbstractVector
    UtilityMatrix(M::AbstractMatrix, bottom_species=nothing, top_species=nothing)
        if

    end
end
=#
