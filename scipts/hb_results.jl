#=
Created on Thurday 20 Feb 2010
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Plot the results of the honeybee spillover.
=#

using StatsPlots, HypothesisTests, BSON

BSON.@load "results/spillover.bson" Ps_during Ps_after M M_a M_b M_ab kl_ot_ab  kl_ot_a kl_ot_b kl_ot_softmax kl_ot_neutral

# tests
# -----

# test if kl_ot_ab < kl_ot_neutral
@show wsr = SignedRankTest(kl_ot_neutral, kl_ot_ab)
@show pvalue(wsr, tail=:right)


# plotting
# --------

p_mut = boxplot(kl_ot_neutral, label="neutral", legend=:none)
boxplot!(kl_ot_ab, label="OT (both)")
boxplot!(kl_ot_a, label="OT (plants fixed)")
boxplot!(kl_ot_b, label="OT (poll. fixed)")
boxplot!(kl_ot_softmax, label="OT (none)")
ylabel!("Kullback-Leibler (nats)")
xticks!(1:5, ["neutral", "OT (both)", "OT (plants fixed)", "OT (poll. fixed)", "OT (none)"])
title!("Predicted versus observed pollination\nduring honeybee spillover")

savefig(p_mut, "figures/honeybee_spillover.svg")

# Same figure, but on training
using OTSIN

kl_ot_ab = Float64[]
kl_ot_a = Float64[]
kl_ot_b = Float64[]
kl_ot_softmax = Float64[]
kl_ot_neutral = Float64[]

for P in Ps_after
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

p_mut = boxplot(kl_ot_neutral, label="neutral", legend=:none)
boxplot!(kl_ot_ab, label="OT (both)")
boxplot!(kl_ot_a, label="OT (plants fixed)")
boxplot!(kl_ot_b, label="OT (poll. fixed)")
boxplot!(kl_ot_softmax, label="OT (none)")
ylabel!("Kullback-Leibler (nats)")
xticks!(1:5, ["neutral", "OT (both)", "OT (plants fixed)", "OT (poll. fixed)", "OT (none)"])
title!("Predicted versus observed pollination\nafter honeybee spillover")

savefig(p_mut, "figures/honeybee_spillover_train.svg")
