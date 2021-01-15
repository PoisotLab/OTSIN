### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ cbdf4e42-2985-11eb-2041-dbd79e434686
using Plots, OTSIN, StatsBase, NamedArrays

# ╔═╡ d94c4784-29bb-11eb-1050-517f3dd80fa0
md"""
Illustration of the optimal transportation on a seed dispersal dataset `M_SD_004`.

```
	
Spindalis portoricensis
Euphonia musica
Nesospingus sp1 M_SD_004
Loxigilla portoricensis
Vireo altiloquous
Melanerpes portoricensis
Columba squamosa
Margarops fuscatus
Tyrannus dominicensis
Zenaida asiatica
Coereba flaveola
Dendroica tigrina
Vireo flavirostris
Dendroica caerulescens
Icterus dominicensis
Myarchus antillarum
Tiaris bicolor
Todus mexicanus
Tyrannus caudifasciatus
Vireo latimeri
Phoradendron sp1 M_SD_004	10	55	3	9	1	0	0	0	0	0	1	0	1	0	0	0	0	0	0	0
Schefflera morototoni	56	1	4	1	3	0	5	1	4	0	1	1	0	0	0	0	0	0	1	0
Cecropia schreberiana	29	0	6	6	4	3	0	1	0	0	2	3	0	1	0	0	0	0	0	0
Ficus sp1 M_SD_004	34	0	4	0	1	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
Guarea guidonia	0	0	5	17	14	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0
Cordia sulcata	28	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Alchornea latifolia	3	0	0	14	3	0	0	0	0	5	0	0	0	0	0	0	0	1	0	1
Miconia serrulata	6	2	12	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Clusia rosea	0	0	1	0	11	3	0	0	1	0	0	0	0	0	1	1	0	0	0	0
Anthurium scandens	1	16	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Miconia affinis	2	0	10	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Inga vera	0	0	8	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Musa acuminata	1	0	3	4	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Inga laurina	4	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Dendropemon bicolor	0	0	3	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Myrcine coriacea	0	0	2	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Syzigium jambos	0	0	1	3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Andira inermis	1	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Coffea arabica	0	0	0	1	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Miconia racemosa	0	0	2	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Palicourea guianensis	0	0	1	2	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Tilandsia sp1 M_SD_004	0	0	3	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Buchenavia capitata	0	0	0	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Casearia arborea	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Citrus sinensis	1	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Dendropanax arboreus	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Henriettea fascicularis	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Momordica charantia	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Palicourea crocea	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Philodendron angustatum	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Piper sp1 M_SD_004	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
Solanum rugosum	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
Unidentified sp1 M_SD_004	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
Zanthoxylum martinicensis	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
```

"""

# ╔═╡ 8dae9448-29bc-11eb-18b0-fb55e9fb1be3
plant_species = ["Phoradendron sp1 M_SD_004",
	"Schefflera morototoni",
	"Cecropia schreberiana",
	"Ficus sp1 M_SD_004",
	"Guarea guidonia",
	"Cordia sulcata",
	"Alchornea latifolia",
	"Miconia serrulata",
	"Clusia rosea",
	"Anthurium scandens",
	"Miconia affinis",
"Inga vera",
"Musa acuminata",
"Inga laurina",
"Dendropemon bicolor",
"Myrcine coriacea",
"Syzigium jambos",
"Andira inermis",
"Coffea arabica",
"Miconia racemosa",
"Palicourea guianensis",
"Tilandsia sp1 M_SD_004",
"Buchenavia capitata",
"Casearia arborea",
"Citrus sinensis",
"Dendropanax arboreus",
"Henriettea fascicularis",
"Momordica charantia",
"Palicourea crocea",
"Philodendron angustatum",
"Piper sp1 M_SD_004",
"Solanum rugosum",
"Unidentified sp1 M_SD_004",
"Zanthoxylum martinicensis"]

# ╔═╡ 67d66742-2b32-11eb-1a58-370b8ca72b53
length(plant_species)

# ╔═╡ 37e622b0-29bc-11eb-2ede-a9430c312543
bird_species = ["Spindalis portoricensis",
	"Euphonia musica",
	"Nesospingus sp1 M_SD_004",
	"Loxigilla portoricensis",
	"Vireo altiloquous",
	"Melanerpes portoricensis",
	"Columba squamosa",
	"Margarops fuscatus",
	"Tyrannus dominicensis",
	"Zenaida asiatica",
	"Coereba flaveola",
	"Dendroica tigrina",
	"Vireo flavirostris",
	"Dendroica caerulescens",
	"Icterus dominicensis",
	"Myarchus antillarum",
	"Tiaris bicolor",
	"Todus mexicanus",
	"Tyrannus caudifasciatus",
	"Vireo latimeri"]

# ╔═╡ 60219008-2b32-11eb-3279-37c1b9489162
length(bird_species)

# ╔═╡ ed1c921c-2986-11eb-28ad-e1f7d3489193
Y = [10 55 3 9 1 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0;
56 1 4 1 3 0 5 1 4 0 1 1 0 0 0 0 0 0 1 0;
29 0 6 6 4 3 0 1 0 0 2 3 0 1 0 0 0 0 0 0;
34 0 4 0 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 5 17 14 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0;
28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
3 0 0 14 3 0 0 0 0 5 0 0 0 0 0 0 0 1 0 1;
6 2 12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 11 3 0 0 1 0 0 0 0 0 1 1 0 0 0 0;
1 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
2 0 10 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 8 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 3 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
4 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 2 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

# ╔═╡ 06c42ff4-29be-11eb-23d3-bffa5ba9a554
n_interactions = sum(Y)

# ╔═╡ 08688ab4-298a-11eb-0ee0-0f00e8ee6248
n, m = size(Y)

# ╔═╡ 3ec1637c-2987-11eb-27c3-7d291ca7306e
a, b = marginals(Y)

# ╔═╡ 4859100e-2987-11eb-3f77-b5825e36a555
P = normalize(Y)

# ╔═╡ 398ed8c4-29bb-11eb-3089-afafbc8345af
plot_Y = marginalsheatmap(Y, cmap=cgrad([:white, :blue]), showaxis=true, xticks=(5:5:m), yticks=5:5:n, xlabel="bird species", ylabel="plant species")

# ╔═╡ 50b53a72-2987-11eb-27b5-a53dc0b26e5d
M = fitM(P, fix_a=true, fix_b=true, γ=1e-4)

# ╔═╡ 88a8f162-2987-11eb-2ed8-d5689f6e2895
plot_M = heatmap(M, color=cgrad([:blue,:white,:red]), xlabel="bird species", ylabel="plant species", title="fitted M")

# ╔═╡ 66caee3c-2988-11eb-188a-ef54b7e65c39
Q = optimaltransport(M, a, b)

# ╔═╡ c72c0ccc-29bb-11eb-1819-f361cf9caa6f
Qpu = optimaltransport(M, ones(n)/n, b)

# ╔═╡ 8de8b87c-29be-11eb-3022-83ee7e3f562b
Qbo = optimaltransport(M, a)

# ╔═╡ 368d73e2-2989-11eb-3dc5-2f89432bff64
plot_Qpu = marginalsheatmap(Qpu, cmap=cgrad([:white, :purple]), showaxis=true, xticks=(5:5:m), yticks=5:5:n, xlabel="bird species", ylabel="plant species")

# ╔═╡ 98e014b6-29be-11eb-0713-1de8527dc3d2
plot_Qbo = marginalsheatmap(Qbo, cmap=cgrad([:white, :purple]), showaxis=true, xticks=(5:5:m), yticks=5:5:n, xlabel="bird species", ylabel="plant species")

# ╔═╡ 265d6474-29bf-11eb-2e20-a1180d19e9e8
p_total = plot(plot_Y, plot_M, plot_Qpu, plot_Qbo, size=(1200, 1000))

# ╔═╡ 7597e448-29c3-11eb-2262-951017f11586
savefig(p_total, "../figures/SD_sim.svg")

# ╔═╡ 7762f2c6-2988-11eb-33b1-a3817c1041b7
mean(abs.(Y .- Q)) 

# ╔═╡ 26ca4e9e-2989-11eb-06ae-fb1ec5ab97a8
kldivergence(P, normalize(Q))

# ╔═╡ Cell order:
# ╟─d94c4784-29bb-11eb-1050-517f3dd80fa0
# ╠═8dae9448-29bc-11eb-18b0-fb55e9fb1be3
# ╠═67d66742-2b32-11eb-1a58-370b8ca72b53
# ╠═37e622b0-29bc-11eb-2ede-a9430c312543
# ╠═60219008-2b32-11eb-3279-37c1b9489162
# ╠═cbdf4e42-2985-11eb-2041-dbd79e434686
# ╠═ed1c921c-2986-11eb-28ad-e1f7d3489193
# ╠═06c42ff4-29be-11eb-23d3-bffa5ba9a554
# ╠═08688ab4-298a-11eb-0ee0-0f00e8ee6248
# ╠═3ec1637c-2987-11eb-27c3-7d291ca7306e
# ╠═4859100e-2987-11eb-3f77-b5825e36a555
# ╠═398ed8c4-29bb-11eb-3089-afafbc8345af
# ╠═50b53a72-2987-11eb-27b5-a53dc0b26e5d
# ╠═88a8f162-2987-11eb-2ed8-d5689f6e2895
# ╠═66caee3c-2988-11eb-188a-ef54b7e65c39
# ╠═c72c0ccc-29bb-11eb-1819-f361cf9caa6f
# ╠═8de8b87c-29be-11eb-3022-83ee7e3f562b
# ╠═368d73e2-2989-11eb-3dc5-2f89432bff64
# ╠═98e014b6-29be-11eb-0713-1de8527dc3d2
# ╠═265d6474-29bf-11eb-2e20-a1180d19e9e8
# ╠═7597e448-29c3-11eb-2262-951017f11586
# ╠═7762f2c6-2988-11eb-33b1-a3817c1041b7
# ╠═26ca4e9e-2989-11eb-06ae-fb1ec5ab97a8
