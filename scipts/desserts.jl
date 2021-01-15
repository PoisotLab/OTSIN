#=
Created on Wednesday 08 July 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Example of dessert distribution to illustrate the concept.
=#

using OTSIN, Plots, StatsPlots, Printf
using LaTeXStrings

a = [3, 3, 3, 2, 2, 2]
b = [5, 4, 3, 3]

@assert sum(a) ≈ sum(b)


desserts = [:merveilleux, :eclair, Symbol("chocolate mousse"), Symbol("carrot cake")]
persons = [:Bernard, :Jan, :Willem, :Laura, :Laure, :Margot]

M = [2 2 1 0;
     0 -2 -2 2;
     1 2 2 -1;
     -1 2 1 -2;
     -1 -1 -2 1;
     -1 0 0 2]

Q = optimaltransport(M, a, b, λ=10.0)

λs = [10^i for i in -2:0.05:1]

utilities = Float64[]
entropies = Float64[]
for λ in λs
        Q = optimaltransport(M, a, b, λ=λ)
        Q ./= sum(Q)  # normalize
        push!(utilities, utility(M, Q))
        push!(entropies, entropy(Q))
end

groupedbar(String.(persons), Q, labels=reshape(desserts, 1, :), bar_position = :stack,
        title="Neutrally-driven")

anim = @animate for λ in λs
  Q = optimaltransport(M, a, b, λ=λ)
  pQ = groupedbar(String.(persons), Q, labels=reshape(desserts, 1, :),
                bar_position = :stack,
                title=@sprintf("Optimal coupling\n lambda=%.2f", λ), ylabel="portions")
  put = plot(λs, utilities, label="utility", lw=2)
  plot!(λs, entropies, label="entropy", lw=2)
  xaxis!(:log)
  xlabel!(L"$\lambda$")
  vline!([λ], label="", lw=2, color=:red, ls=:dash)
  plot(pQ, put, size=(800, 500))
end

gif(anim, "figures/desserts.gif", fps = 5)
