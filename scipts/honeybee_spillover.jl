#=
Created on Thursday 2 April 2020
Last update:

@author: Michiel Stock
michielfmstock@gmail.com

Experiments on the honeybee spillover dataset.
=#

using DataFrames, CSV
using BSON
using Plots, Statistics

interactions = CSV.read("data/Magrach2017_network.csv")
first(interactions, 5)

names(interactions)

sites = unique(interactions.site_id)

pollinators = unique(interactions.pollinator_species);
plants = unique(interactions.plant_species);

make_inverse_index = l -> Dict(n => i for (i, n) in enumerate(l))

pollinators_inv_index = make_inverse_index(pollinators);
plants_inv_index = make_inverse_index(plants);

println("Plants:")
for p in plants
    println("- $p")
end

println("Pollinators:")
for p in pollinators
    println("- $p")
end

function make_incidence(interactions, pl_inv_index, pol_inv_index)
    Y = zeros(length(pl_inv_index), length(pol_inv_index))
    for (pl, pol, int) in zip(interactions[:,:plant_species], interactions[:,:pollinator_species], interactions[:,:transect])
        i, j = pl_inv_index[pl], pol_inv_index[pol]
        Y[i,j] += int
    end
    return Y
end

Y = make_incidence(interactions, plants_inv_index, pollinators_inv_index)

networks_after = Dict{String, Array{Float64,2}}()
networks_during = Dict{String, Array{Float64,2}}()


for site in interactions.site_id
    interactions_site = interactions[interactions.site_id.==site,:]
    Y_after = make_incidence(interactions_site[interactions_site.period.=="after",:],
                            plants_inv_index, pollinators_inv_index)
    Y_during = make_incidence(interactions_site[interactions_site.period.=="during",:],
                            plants_inv_index, pollinators_inv_index)

    networks_after[site] = Y_after
    networks_during[site] = Y_during
end

for (time, networks) in zip(["after", "during"], [networks_after, networks_during])
    println("$time")
    for (site, Y) in networks
        println("\t- $site: $(sum(Y)) ($(sum(Y.>0)))")
    end
end

using OTSIN

Ps = [(Y .+ 1e-4) |> OTSIN.normalize for Y in values(networks_after)];

γ = 0.010

@time M_ab = fitM(Ps, fix_a=true, fix_b=true, maxitter=20, γ=γ);

@time M_a = fitM(Ps, fix_a=true, fix_b=false, γ=γ);

@time M_b = fitM(Ps, fix_a=false, fix_b=true, γ=γ);

@time M = fitM(Ps, fix_a=false, fix_b=false, γ=γ);

plot(heatmap(M, title="M"), heatmap(M_a, title="M (a fixed)"),
        heatmap(M_b, title="M (b fixed)"), heatmap(M_ab, title="M (a, b fixed)"))
savefig("figures/M_honeybee_spillover.svg")

Ps_during = [Y |> OTSIN.normalize for Y in values(networks_during)];

kl_ot_ab = Float64[]
kl_ot_a = Float64[]
kl_ot_b = Float64[]
kl_ot_softmax = Float64[]
kl_ot_neutral = Float64[]

for P in Ps_during
    a, b = marginals(P)
    P = P[a.>0, b.>0]
    M = M_ab[a.>0, b.>0]
    a, b = marginals(P)
    # optimal transport both
    Q = optimaltransport(M, a, b)
    push!(kl_ot_ab, KL(P, Q))
    Q .= optimaltransport(M, a, nothing)
    push!(kl_ot_a, KL(P, Q))
    Q .= optimaltransport(M, nothing, b)
    push!(kl_ot_b, KL(P, Q))
    Q .= optimaltransport(M, nothing, nothing)
    push!(kl_ot_softmax, KL(P, Q))
    Q .= optimaltransport(a, b)
    push!(kl_ot_neutral, KL(P, Q))
end

Ps_after = Ps

BSON.@save "results/spillover.bson" Ps_during Ps_after M M_a M_b M_ab kl_ot_ab  kl_ot_a kl_ot_b kl_ot_softmax kl_ot_neutral γ
