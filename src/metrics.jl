"""
    entropy(P)

Computes the entropy of a distribution.
"""
function entropy(P; norm=false)
    if norm
        P /= sum(P)
    end
    return - sum(P[P.>0] .* log.(P[P.>0]))
end

"""
    relative_entropy(P, Q)

Computes the relative entropy of P to Q.
"""
relative_entropy(P, Q) = - sum(P .* log.(Q))

"""
    KL(P, Q)

Computes the Kullback-Leibler divergence of P to Q
"""
KL(P, Q) = relative_entropy(P, Q) - entropy(P)


"""
    utility(P, M)

Computes the average utility.
"""
function utility(P, M)
    return sum(P .* M)
end

function report_utility(P, M, bottom_species=nothing, top_species=nothing)
    @assert size(P) == size(M) "P and M should be of the same size"
    n, m = size(P)
    bottom_species == nothing && bottom_species == 1:n
    top_species == nothing && top_species == 1:m
    U = P .* M
    a, b = marginals(P)
    println("Bottom species")
    for i in 1:n
        println("\t$(bottom_species[i]): utility = $(sum(U[i,:]) / a[i]), entropy = $(
                                entropy(P[i,:] |> x -> x / sum(x)))")
    end
    println()
    println("Top species")
    for j in 1:m
        println("\t$(top_species[j]): utility = $(sum(U[:,j] / b[j])), entropy = $(
                                entropy(P[:,j] |> x -> x / sum(x)))")
    end
    println()
    println("Average utility = $(sum(U))")
    println("Average entropy = $(entropy(P))")
end
