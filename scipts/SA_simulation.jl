### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 04911b34-229d-11eb-1658-3d8632a4e006
using Distributions, StatsBase, Plots, OTSIN

# ╔═╡ bc3b1416-2355-11eb-16b4-c1e191b1dc94
md"""
# OT solutions arise from local choices

In this notebook, I will show that the optimal transportation model arises from local choises of species interactions. It suffices to do local choices of species based on the **utility** and Metropolis heuristic. This leads to behaviour for which the optimal transportation problem is an exact solution.
"""

# ╔═╡ 8547775c-2359-11eb-1a41-673a988bfdca
md"First, we define a utility matrix."

# ╔═╡ 6144ef0e-2333-11eb-17a9-0392d858553e
begin
	Ula = zeros(30, 25)
	Ula[1:10, 1:12] .= 2.0
	Ula[9:20, 10:20] .= 1.5
	Ula[22:30, 20:end] .= 1.0
	Ula[11:end, 1:10] .= 0.5
	heatmap(Ula)
end

# ╔═╡ d4fe37fa-2335-11eb-31f6-bfda1afcd852
Usm = [3 3 3 1 0;
	   3 2.5 2 1 0.5;
	   2.5 2 2 0 0;
       1  1 0 0 1;
       1 0 0 0.5 1]

# ╔═╡ 60672f2a-229d-11eb-1f7f-a59e191a8a65
#M = 10*[ones(15, 10) zeros(15, 10); 0.1ones(10, 12) 0.5ones(10, 8)];
M = Usm

# ╔═╡ 3a16e734-229d-11eb-3187-d721246abe92
n, m = size(M)

# ╔═╡ bd664532-229d-11eb-19aa-1d5847ad8687
a = ones(n) / n

# ╔═╡ c99500bc-229d-11eb-1459-fb202f664fba
b = range(10, 1, length=m) |> x -> x/sum(x)

# ╔═╡ eac665fa-229d-11eb-1404-73687f0ea99e
multinomr(n, p) = rand(Multinomial(n, p))

# ╔═╡ 51c9802a-229e-11eb-0800-d180f7c63640
sample(1:n, Weights(a), 2)

# ╔═╡ 214d7c6c-229e-11eb-2b54-cf6068a5ff69

function simulated_annealing!(P, M, λ; nsteps=10000)
	n, m = size(P)
	n_interactions = sum(P)
	rowweigths = Weights(sum(P, dims=2)[:])
	for t in 1:nsteps
		i1, i2 = sample(1:n, rowweigths, 2)
		j1 = sample(1:m, Weights(P[i1,:]))
		j2 = sample(1:m, Weights(P[i2,:]))
		ΔM = (M[i1,j2] + M[i2, j1]) - (M[i1,j1] + M[i2,j2])
		if ΔM > 0.0 || rand() ≤ exp(ΔM * λ)
			P[i1,j1] -= 1
			P[i2,j2] -= 1
			P[i1,j2] += 1
			P[i2,j1] += 1
		end
	end
	return P
end


# ╔═╡ b85f7bca-229f-11eb-1e15-c7cb10437450
n_interactions = 200

# ╔═╡ c0c50460-229f-11eb-3a58-7d43a336d5e3
a_sampled = multinomr(n_interactions, a)

# ╔═╡ ce10d4b6-229f-11eb-06bd-3720e2a99287
b_sampled = multinomr(n_interactions, b)

# ╔═╡ d9bc7310-229f-11eb-37fe-3f30f4814257
function findP(a_sampled, b_sampled, n_interactions)
	n, m = length(a_sampled), length(b_sampled)
	P = zeros(Int, n, m)
	while n_interactions > 0
		i, j = rand(1:n), rand(1:m)
		if sum(P[i,:]) < a_sampled[i] && sum(P[:,j]) < b_sampled[j]
			P[i,j] += 1
			n_interactions -= 1
		end
	end
	return P
end

# ╔═╡ 7fb496e4-22a0-11eb-1080-ebf1eb4f42bc
P = findP(a_sampled, b_sampled, n_interactions)

# ╔═╡ fd192064-22a0-11eb-1fcc-b11f49ad44ab
simulated_annealing!(P, M, 10; nsteps=1000_000)

# ╔═╡ 2f8e074e-22a1-11eb-136d-2b83877dd889
marginalsheatmap(P)

# ╔═╡ 6896da6a-22a2-11eb-1fcc-839ef4100dc1
Pot = optimaltransport(M, a, b, λ=10, maxitter=100)

# ╔═╡ 806ed552-22a2-11eb-203a-419eee9307a3
marginalsheatmap(Pot)

# ╔═╡ 295134bc-22a3-11eb-3e28-fb55a02936e8
λs = (10).^(-2:0.2:1.5)

# ╔═╡ 76315a66-22ab-11eb-0a31-adee0035db0d
n_experiments = length(λs)

# ╔═╡ 844fa4c2-22ab-11eb-1186-f1b3644cc72f
n_repetitions = 50

# ╔═╡ de07a940-22ab-11eb-0113-83f986874b02
numb_interactions = [100, 1000, 10_000]

# ╔═╡ bf4dc6c6-236e-11eb-17b5-63f843e8cce7
minheatmap(P; kwargs...) = heatmap(P, legend=false, seriescolor=cgrad([:white, :blue]),showaxis=true, xticks=(1:size(P,1), ('a':'z')[1:size(P,1)]),
				yticks=(1:size(P,2), ('A':'Z')[1:size(P,2)]); kwargs...)

# ╔═╡ b6270016-229d-11eb-1d5a-1b124847f727
pM = minheatmap(M, color=cgrad([:white, :red]))

# ╔═╡ f8797b2a-27ea-11eb-2be7-d76ac2a24e4c
savefig(pM, "../figures/M_sa_sim.svg")

# ╔═╡ c1c25be2-236e-11eb-0b1e-27ce2c4a7416
minheatmap(P)

# ╔═╡ 7fba27a0-22a5-11eb-1a5f-85c0b7fe545c
function entropy(P)
	total = sum(P)
	logtotal = log(total)
	h = 0.0
	for Pij in P
		if Pij > 0
			Pij /= total
			h += -Pij * log(Pij)
		end
	end
	return h
end
			

# ╔═╡ f51dd178-22a7-11eb-286d-756e77b64941
average_utity(P, M) = sum(P .* M) / sum(P)

# ╔═╡ e29aaa28-22b4-11eb-097b-2379ac1d3ee8
average_utity(P, M)

# ╔═╡ f4b75eca-22a2-11eb-3454-b7cf0c58ae0b
begin
	av_utilities = zeros(n_experiments, length(numb_interactions));
for (ind_int, n_int) in enumerate(numb_interactions)
	for rep in 1:n_repetitions
		P = findP(multinomr(n_int, a), multinomr(n_int, b), n_int)
		for (i, λ) in enumerate(λs)
			simulated_annealing!(P, M, λ; nsteps=100000)
			av_utilities[i, ind_int] += average_utity(P, M)
		end
		
	end
		
end
av_utilities ./= n_repetitions
end

# ╔═╡ 2dc31e48-2372-11eb-00c0-a3fdd9d08e56
begin
	utilities = zeros(n_experiments, n_repetitions)
	for rep in 1:n_repetitions
		P = findP(multinomr(n_interactions, a), multinomr(n_interactions, b), n_interactions)
		for (i, λ) in enumerate(λs)
			simulated_annealing!(P, M, λ; nsteps=100_000)
			utilities[i, rep] += average_utity(P, M)
		end
	end
	utilities
end

# ╔═╡ 3e1c6016-2374-11eb-3b85-3b419741ec61
μ_ut = mean(utilities, dims=2)

# ╔═╡ 459210c0-2374-11eb-1fab-a123f5b4fe1f
σ_ut = std(utilities, dims=2)

# ╔═╡ 73810eaa-22af-11eb-2417-6d34b137efac
ot_utilites_theo = [average_utity(optimaltransport(M, a, b, λ=λ), M) for λ in λs]

# ╔═╡ 846a5720-22ad-11eb-1adf-d5f29a11edd2
begin
	plambda = plot(λs, [μ_ut μ_ut], fillrange=[μ_ut+1σ_ut μ_ut-1σ_ut], fillalpha=0.4, color=:orange, label=["sampled utitily" ""], lw=4, xscale=:log, xlabel="lambda", ylabel="average utility")
	#plot!(λs, μ_ut, label = "", , color=:orange)
	plot!(λs, ot_utilites_theo, label="OT bound", ls=:dash, lw=4, color=:black)
	plot!(legend=:topleft)
end

# ╔═╡ e91cf2a0-236c-11eb-13b7-99f464bd8c2c
begin
	l = @layout [a{0.60h}; b c d; e f g]
	p_sum = plot(plambda,
		minheatmap(optimaltransport(M, a, b, λ=0.05), title="lambda=0.05"),
		minheatmap(optimaltransport(M, a, b, λ=0.5), title="lambda=0.5"),
		minheatmap(optimaltransport(M, a, b, λ=5), title="lambda=5"),
		minheatmap(simulated_annealing!(P, M, 0.05, nsteps=1000_000), seriescolor=cgrad([:white, :green])),
		minheatmap(simulated_annealing!(P, M, 0.5, nsteps=1000_000), seriescolor=cgrad([:white, :green])),
		minheatmap(simulated_annealing!(P, M, 5, nsteps=1000_000), seriescolor=cgrad([:white, :green])),
		layout=l, size=(600, 800))
	savefig("../figures/OT_MH_approx.svg")
	p_sum
end

# ╔═╡ 594a8728-2374-11eb-01bc-e9bf600f89d3
begin
	pband = plot(λs, av_utilities,
			labels = reshape(["n interactions=$nintr" for nintr in numb_interactions], 1, :), lw=2, xscale=:log, xlabel="lambda", ylabel="average utility")
	plot!(λs, ot_utilites_theo, label="OT bound", ls=:dash, lw=4, color=:black)
	plot!(legend=:topleft)
end

# ╔═╡ 2184453a-22a8-11eb-301f-9553f16a0c4c
sum(Pot)

# ╔═╡ ac9d6a72-22b5-11eb-2060-21112685115f
Ps = reshape(rand(Multinomial(10, vec(Pot))), n, m)

# ╔═╡ 28543ba0-22b6-11eb-3457-937945c42c47
average_utity(Ps, M)

# ╔═╡ Cell order:
# ╠═bc3b1416-2355-11eb-16b4-c1e191b1dc94
# ╠═04911b34-229d-11eb-1658-3d8632a4e006
# ╟─8547775c-2359-11eb-1a41-673a988bfdca
# ╠═3a16e734-229d-11eb-3187-d721246abe92
# ╠═6144ef0e-2333-11eb-17a9-0392d858553e
# ╠═b6270016-229d-11eb-1d5a-1b124847f727
# ╠═f8797b2a-27ea-11eb-2be7-d76ac2a24e4c
# ╠═d4fe37fa-2335-11eb-31f6-bfda1afcd852
# ╠═60672f2a-229d-11eb-1f7f-a59e191a8a65
# ╠═bd664532-229d-11eb-19aa-1d5847ad8687
# ╠═c99500bc-229d-11eb-1459-fb202f664fba
# ╠═eac665fa-229d-11eb-1404-73687f0ea99e
# ╠═51c9802a-229e-11eb-0800-d180f7c63640
# ╠═214d7c6c-229e-11eb-2b54-cf6068a5ff69
# ╠═b85f7bca-229f-11eb-1e15-c7cb10437450
# ╠═c0c50460-229f-11eb-3a58-7d43a336d5e3
# ╠═ce10d4b6-229f-11eb-06bd-3720e2a99287
# ╠═d9bc7310-229f-11eb-37fe-3f30f4814257
# ╠═7fb496e4-22a0-11eb-1080-ebf1eb4f42bc
# ╠═fd192064-22a0-11eb-1fcc-b11f49ad44ab
# ╠═2f8e074e-22a1-11eb-136d-2b83877dd889
# ╠═e29aaa28-22b4-11eb-097b-2379ac1d3ee8
# ╠═6896da6a-22a2-11eb-1fcc-839ef4100dc1
# ╠═806ed552-22a2-11eb-203a-419eee9307a3
# ╠═295134bc-22a3-11eb-3e28-fb55a02936e8
# ╠═76315a66-22ab-11eb-0a31-adee0035db0d
# ╠═844fa4c2-22ab-11eb-1186-f1b3644cc72f
# ╠═de07a940-22ab-11eb-0113-83f986874b02
# ╠═f4b75eca-22a2-11eb-3454-b7cf0c58ae0b
# ╠═2dc31e48-2372-11eb-00c0-a3fdd9d08e56
# ╠═3e1c6016-2374-11eb-3b85-3b419741ec61
# ╠═459210c0-2374-11eb-1fab-a123f5b4fe1f
# ╠═73810eaa-22af-11eb-2417-6d34b137efac
# ╠═846a5720-22ad-11eb-1adf-d5f29a11edd2
# ╠═594a8728-2374-11eb-01bc-e9bf600f89d3
# ╠═e91cf2a0-236c-11eb-13b7-99f464bd8c2c
# ╠═bf4dc6c6-236e-11eb-17b5-63f843e8cce7
# ╠═c1c25be2-236e-11eb-0b1e-27ce2c4a7416
# ╠═7fba27a0-22a5-11eb-1a5f-85c0b7fe545c
# ╠═f51dd178-22a7-11eb-286d-756e77b64941
# ╠═2184453a-22a8-11eb-301f-9553f16a0c4c
# ╠═ac9d6a72-22b5-11eb-2060-21112685115f
# ╠═28543ba0-22b6-11eb-3457-937945c42c47
