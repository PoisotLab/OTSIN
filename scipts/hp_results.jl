#=
Created on Thurday 20 Feb 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Process the results of the host-parasite data.
=#

using StatsPlots, HypothesisTests, BSON, NamedArrays

BSON.@load "results/hp.bson" Pstrain Pstest M kl_ot_ab  kl_ot_a kl_ot_b kl_ot_softmax kl_ot_neutral γ

# HYPOTHESIS TESTS
# ----------------

Δkl = kl_ot_neutral - kl_ot_ab

mean(Δkl) / (std(Δkl) / √(length(Δkl)))

@show wsr = SignedRankTest(kl_ot_neutral, kl_ot_ab)
@show pvalue(wsr, tail=:right)

# FIGURES
# -------


p_hp = boxplot(kl_ot_neutral, label="neutral", legend=:none)
boxplot!(kl_ot_ab, label="OT (both)")
boxplot!(kl_ot_a, label="OT (parasites fixed)")
boxplot!(kl_ot_b, label="OT (hosts fixed)")
boxplot!(kl_ot_softmax, label="OT (none)")
ylabel!("Kullback-Leibler (nats)")
xticks!(1:5, ["neutral", "OT (both)", "OT (parasites fixed)", "OT (hosts fixed)", "OT (none)"])
title!("Predicted versus observed coupling\nhost-parasite networks")

savefig(p_hp, "figures/hostparasite.svg")

# on the training data

kl_ot_ab = Float64[]
kl_ot_a = Float64[]
kl_ot_b = Float64[]
kl_ot_softmax = Float64[]
kl_ot_neutral = Float64[]

cutoff = 1e-6  # remove sp. with no interactions

for P in Pstrain
    P = Matrix(P)
    P ./=sum(P)
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

p_hp = boxplot(kl_ot_neutral, label="neutral", legend=:none)
boxplot!(kl_ot_ab, label="OT (both)")
boxplot!(kl_ot_a, label="OT (parasites fixed)")
boxplot!(kl_ot_b, label="OT (hosts fixed)")
boxplot!(kl_ot_softmax, label="OT (none)")
ylabel!("Kullback-Leibler (nats)")
xticks!(1:5, ["neutral", "OT (both)", "OT (parasites fixed)", "OT (hosts fixed)", "OT (none)"])
title!("Predicted versus observed coupling\nhost-parasite networks (training)")

savefig(p_hp, "figures/hostparasite_train.svg")
