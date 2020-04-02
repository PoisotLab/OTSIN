#=
Created on Thursday 2 April 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Experiments on the host-parasite datasets.
=#

using DataFrames, CSV
using Plots, Statistics
using NamedArrays

using EcologicalNetworks

all_hp_data = filter(x -> occursin("Hadfield", x.Reference), web_of_life());

ids = getfield.(all_hp_data, :ID);
networks = web_of_life.(ids)

N = reduce(union, convert.(BinaryNetwork, networks))

parasites, hosts = species(N, dims=1), species(N, dims=2);

M = NamedArray(zeros(size(N)...), (parasites, hosts))

networks_completed = NamedArray{Float64,2}[]

for (network, id) in zip(networks, ids)
    Y = similar(M)
    fill!(Y, 0.0)

    for p in species(network, dims=1)
        for h in species(network, dims=2)
            Y[p, h] = network[p, h]
        end
    end
    push!(networks_completed, Y)
end

networks_completed

sum.(networks_completed)

using OTSIN, StatsBase

train_inds = sample(1:51, 26, replace=false)
test_inds = setdiff(1:51, train_inds);

γ = 1e-2

Ps = [i in train_inds ? Y .+ 1e-4 : Y for (i, Y) in enumerate(networks_completed)]
Ps = OTSIN.normalize.(Ps)
@time M = fitM(Matrix.(Ps[train_inds]), fix_a=true, fix_b=true, maxitter=5, γ=γ)
#@time M = fitM(Matrix.(Ps[train_inds]), fix_a=false, fix_b=false, maxitter=10, γ=1.0)

heatmap(M, color=:pu_or)
title!("Global utility matrix host-parasite")
xlabel!("hosts")
ylabel!("parasites")

savefig("figures/hostparM.svg")

kl_ot_ab = Float64[]
kl_ot_a = Float64[]
kl_ot_b = Float64[]
kl_ot_softmax = Float64[]
kl_ot_neutral = Float64[]

cutoff = 1e-6  # remove sp. with no interactions

for P in Ps[test_inds]
    P = Matrix(P)
    a, b = marginals(P)
    # remove empty rows/columns
    P = P[a.>cutoff, b.>cutoff]
    Mss = M[a.>cutoff, b.>cutoff]
    a, b = marginals(P)
    # optimal transport both
    Q = optimaltransport(Mss, a, b)
    push!(kl_ot_ab, KL(P, Q))
    Q .= optimaltransport(Mss, a, nothing)
    push!(kl_ot_a, KL(P, Q))
    Q .= optimaltransport(Mss, nothing, b)
    push!(kl_ot_b, KL(P, Q))
    Q .= optimaltransport(Mss, nothing, nothing)
    push!(kl_ot_softmax, KL(P, Q))
    Q .= optimaltransport(a, b)
    push!(kl_ot_neutral, KL(P, Q))
end

using BSON

Ps = convert.(Array, Ps)
Pstrain = Ps[train_inds]
Pstest = Ps[test_inds]

BSON.@save "results/hp.bson" Pstrain Pstest M kl_ot_ab  kl_ot_a kl_ot_b kl_ot_softmax kl_ot_neutral γ
