#=
Created on Thursday 2 April 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Simple simulation to illustrate OTSIN.
=#


using OTSIN
using Plots
using StatsPlots
using Distributions
using DataFrames

# GENERATE DATA
# -------------

top_traits = [Normal(1.6, 0.4), Normal(2.2, 0.6), Normal(2.7, 0.3), Normal(5, 0.9)]
top_names = ["A" "B" "C" "D"]

bottom_traits = [1.5, 2, 2.5, 4, 5.4]
bottom_names = ["a", "b", "c", "d", "e"]

gen_inv_trait = Normal(3.7, 1.4)
spec_inv_trait = Normal(2.85, 0.25)

a, b = [0.20, 0.20, 0.1, 0.2, 0.3], [0.30, 0.4, 0.1, 0.2]

@assert sum(a) ≈ sum(b)

# PLOTTING
# --------

xvals = 0.5:0.01:7

top_prefs = [pdf.(top_traits[1], xvals) pdf.(top_traits[2], xvals) pdf.(top_traits[3], xvals) pdf.(top_traits[4], xvals)]

p_traits = plot(xvals, top_prefs, color=:orange, linewidth=2, label="", linestyle=:dash)
vline!(bottom_traits, color=:green, linewidth=2, label="bottom sp.")
vline!(mean.(top_traits), color=:orange, linewidth=2, alpha=0.5, linestyle=:dash, label="top sp.")

# gen invader
plot!(xvals, pdf.(gen_inv_trait, xvals), color=:purple, linewidth=2, label="", linestyle=:dot)
vline!([mean(gen_inv_trait)], color=:purple, linewidth=2, alpha=0.5, linestyle=:dot, label="top gen. invader sp.")

# spec invader
plot!(xvals, pdf.(spec_inv_trait, xvals), color=:red, linewidth=2, label="", linestyle=:dashdot)
vline!([mean(spec_inv_trait)], color=:red, linewidth=2, alpha=0.5, linestyle=:dashdot, label="top spec. invader sp.")


xticks!([bottom_traits; mean.(top_traits); mean.([gen_inv_trait, spec_inv_trait])],
                        [bottom_names; vec(top_names); "E"; "F"])
ylabel!("utility")

p_bar_bottom = bar(bottom_names, a, color=:green, label="", ylabel="frequency", title="Density bottom sp.")
p_bar_top = bar(vec(top_names), b, color=:orange, label="", ylabel="frequency", title="Density top sp.")

# OT SIMULATIONS
# --------------

function make_M(bottom_traits, top_traits; α=10)
    n, m = size(bottom_traits, 1), size(top_traits, 1)
    M = zeros(n, m)

    for i in 1:n
        for j in 1:m
            M[i,j] = pdf(top_traits[j], bottom_traits[i])
        end
    end
    return M
end

M = make_M(bottom_traits, top_traits, α=5)

p_heatmap_M = heatmap(M, xticks=(1:4, top_names), yticks=(1:5, bottom_names), title="Utility matrix M")

l = @layout [a{0.7w}  [b; c; d]]
plot(p_traits, p_heatmap_M, p_bar_top, p_bar_bottom, layout=l, size=(800, 500))
savefig("figures/simulation_outline.svg")

av_utilities = Matrix{Union{Nothing, Float32}}(nothing, 6, 6)

Ql = optimaltransport(M, a, b, λ=1)
Qh = optimaltransport(M, a, b, λ=20)

p_ot_neutral = groupedbar(bottom_names, Ql, labels=top_names, bar_position = :stack,
        title="Neutrally-driven", legend=:none)
p_ot_utility = groupedbar(bottom_names, Qh, labels=top_names, bar_position = :stack,
        title="Utility-driven", legend=:none)

av_utilities[1,1:4] = sum(Qh .* M, dims=1)[:]
av_utilities[2,1:4] = sum(Ql .* M, dims=1)[:]

plot(p_ot_neutral, p_ot_utility)

report_utility(Ql, M, bottom_names, top_names)

report_utility(Qh, M, bottom_names, top_names)

a_removal = copy(a)
a_removal[2] = 0.0
a_removal ./ sum(a_removal)

Q_extinction = optimaltransport(M, a_removal, b, λ=20)

av_utilities[4,1:4] = sum(Q_extinction .* M, dims=1)[:]

p_extincition = groupedbar(bottom_names, Q_extinction, labels=top_names, bar_position = :stack,
        title="Extinction of sp. b", legend=:none)

Q_extinction

M_gen_inv = make_M(bottom_traits, [top_traits; gen_inv_trait])

top_names_inv = [top_names..., "E"]

heatmap(M_gen_inv, xticks=(1:5, top_names_inv), yticks=(1:5, bottom_names))

b_inv = OTSIN.normalize([b; 0.1])

Q_gen_inv = optimaltransport(M_gen_inv, a, b_inv, λ=20.0)

av_utilities[5,1:5] = sum(Q_gen_inv .* M_gen_inv, dims=1)[:]

p_gen_inv = groupedbar(bottom_names, Q_gen_inv, labels=[top_names "E"], bar_position = :stack,
        title="Effect of generalist invader", legend=:top)

report_utility(Q_gen_inv, M_gen_inv, bottom_names, top_names_inv)

M_spec_inv = make_M(bottom_traits, [top_traits; spec_inv_trait])

top_names_inv = [top_names "F"]

heatmap(M_spec_inv, xticks=(1:5, top_names_inv), yticks=(1:5, bottom_names))

Q_spec_inv = optimaltransport(M_spec_inv, a, b_inv, λ=20.0)

av_utilities[6,[1, 2, 3, 4, 6]] = sum(Q_spec_inv .* M_spec_inv, dims=1)[:]

p_spec_inv =  groupedbar(bottom_names, Q_spec_inv, labels=[top_names "F"], bar_position = :stack,
        title="Effect of specialist invader", legend=:top)

plot(p_ot_utility, p_ot_neutral, p_optimal_free, p_extincition, p_gen_inv, p_spec_inv,
        layout=(2, 3), size=(1500, 1.5*500))

savefig("figures/simulation_settings.svg")

p_ut_a, p_ut_b = pert_utility(M, a, b, λ=20)
p_entr_a, p_entr_b = pert_entropy(M, a, b, λ=20)

plot_pert_bottom = scatter(bottom_names, [p_ut_a p_entr_a],
    label=["relative effect utility", "relative effect entropy"],
    title="Perturbations bottom sp.", markershape=[:circle :utriangle])
plot_pert_top = scatter(top_names, [p_ut_b p_entr_b],
    label=["relative effect utility", "relative effect entropy"],
    legend=:none, title="Perturbations top sp.", markershape=[:circle :utriangle])

plot(plot_pert_bottom, plot_pert_top, layout=grid(2, 1), size=(600, 400))
savefig("figures/perturbation.svg")

results = DataFrame(av_utilities, [:A, :B, :C, :D, :E, :F])

println(results)
