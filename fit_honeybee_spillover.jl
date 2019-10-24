using DataFrames, CSV
using Plots

interactions = CSV.read("data/Magrach2017_network.csv")
first(interactions, 5)

sites = unique(interactions.site_id)

pollinators = unique(interactions.pollinator_species);
plants = unique(interactions.plant_species);

make_inverse_index = l -> Dict(n => i for (i, n) in enumerate(l))

pollinators_inv_index = make_inverse_index(pollinators);
plants_inv_index = make_inverse_index(plants);

function make_incidence(interactions, pl_inv_index, pol_inv_index)
    Y = zeros(length(pl_inv_index), length(pol_inv_index))
    for (pl, pol, int) in zip(interactions[:plant_species], interactions[:pollinator_species], interactions[:transect])
        i, j = pl_inv_index[pl], pol_inv_index[pol]
        Y[i,j] += int
    end
    return Y
end

networks_during = Dict{String, Array{Float64,2}}()
networks_after = Dict{String, Array{Float64,2}}()


for site in interactions.site_id
    interactions_site = interactions[interactions.site_id.==site,:]
    Y_during = make_incidence(interactions_site[interactions_site.period.=="during",:],
                            plants_inv_index, pollinators_inv_index)
    Y_after = make_incidence(interactions_site[interactions_site.period.=="after",:],
                            plants_inv_index, pollinators_inv_index)

    networks_during[site] = Y_before
    networks_after[site] = Y_after
end
