#=
Created on Wednesday 9 Janurary 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Trait matching model based on data from Olito et al. 2015.
=#

using CSV, Plots, LaTeXStrings
using OTSIN

interactions = CSV.read("data/Olito2015_SM_A2_rawdata_full/Aggregate network-Table 1.csv", header=2)

plant_morph = CSV.read("data/Olito2015_SM_A2_rawdata_full/Plant morphology-Table 1.csv")
pol_morph = CSV.read("data/Olito2015_SM_A2_rawdata_full/Pollinator morphology-Table 1.csv")

plant_densities = CSV.read("data/Olito2015_SM_A2_rawdata_full/Plant densities-Table 1.csv")


Y = Matrix(interactions[:,2:end])
P = Y / sum(Y)

# match of class
M = Matrix(plant_morph[:,2:5]) * Matrix(pol_morph[:,2:5])'

a, b = marginals(P)

histplants = histogram(sum(Y, dims=2)[:], color=:green, title="interactions plants",
    xlabel="number of interactions", ylabel=:counts, label="");
histpols = histogram(sum(Y, dims=1)[:], color=:yellow,  title="interactions pollinators",
    xlabel="number of interactions", ylabel=:counts, label="");
pmarg = plot(histplants, histpols, layout=(2,1))

λs = 2 .^(-10:0.1:0)

klvals = [KL(P, optimaltransport(M, a, b, λ=λ)) for λ in λs]

λopt = λs[argmin(klvals)]

pλ = plot(λs, klvals, xscale=:log10, lw=2, label="optimal transport")
ylabel!("Kullback-Leibler\ndivergence (nats)")
xlabel!("\$ \\lambda \$")

Q = optimaltransport(M, a, b, λ=λopt)

kl_optimal = minimum(klvals)
kl_neutral = KL(P, optimaltransport(a,b))

hline!([kl_neutral], lw=2, label="neutral", ls=:dash)

pM = heatmap(M, xlabel="plant sp.", ylabel="pol. sp.", title="M", color=:grays) #, colorbar=:none
pP = heatmap(P, xlabel="plant sp.", ylabel="pol. sp.", title="P", color=:heat)
pQ = heatmap(Q, xlabel="plant sp.", ylabel="pol. sp.", title="Q", color=:heat)


l = @layout [  a{0.5w} b{0.5w};
               c{0.3w} d{0.7w}]
plot(pP, pQ, pM, pλ, figsize=(300, 600))

savefig("figures/traitmatch.pdf")
