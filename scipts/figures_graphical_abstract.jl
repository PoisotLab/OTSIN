
ext = ".png"

using OTSIN, LaTeXStrings, Plots

using Colors

lightgrey = RGB(0.98, 0.98, 0.98)

simpleheatmap(M; kwargs...) = heatmap(M, xticks=false,
                yticks=false,colorbar=false; kwargs...)

Y =    [52    7    3    1    1    3    0    1    0    0    0
    18    16    3    10    7    0    0    1    0    0    1
    26    7    13    0    1    3    0    0    0    0    0
    9    11    5    0    0    0    3    0    0    0    0
    13    2    6    2    2    1    0    0    1    1    0
    21    0    1    1    2    0    0    0    0    0    0
    6    3    1    1    0    0    2    0    0    0    0
    0    0    3    0    0    4    0    0    0    0    0
    3    2    0    0    0    0    0    0    0    0    0
    1    2    1    0    0    0    0    0    0    0    0
    0    1    0    2    0    0    0    0    0    0    0
    1    1    0    1    0    0    0    0    0    0    0
    2    0    0    0    0    0    0    0    0    0    0]

P = normalize(Y)

a, b = marginals(P)

plot(bar(a, xticks=false, yticks=false, color="#2a9d8f", label=""),
    bar(b, xticks=false, yticks=false, color="#f4a261", label=""), layout=(2,1))
savefig("figures/GA/marginals$ext")

M = fitM(P)

simpleheatmap(P, color=cgrad([lightgrey, :blue]))
savefig("figures/GA/P$ext")


simpleheatmap(M, color=cgrad([lightgrey, :red]))
savefig("figures/GA/M$ext")

Q = optimaltransport(M, a, b)
simpleheatmap(Q, color=cgrad([lightgrey, :purple]))
savefig("figures/GA/Q$ext")

marginalsheatmap(Q, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qwithmarginals$ext")

Qlow = optimaltransport(M, a, b, λ=0.00001)
marginalsheatmap(Qlow, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qlow$ext")

Qhigh = optimaltransport(M, a, b, λ=100)
marginalsheatmap(Qhigh, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qhigh$ext")

Qa = optimaltransport(M, a, λ=3)
marginalsheatmap(Qa, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qa$ext")

Qb = optimaltransport(M', b, λ=3)'
marginalsheatmap(Qb, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qb$ext")

Qfree = optimaltransport(M, λ=3)
marginalsheatmap(Qfree, cmap=cgrad([lightgrey, :purple]))
savefig("figures/GA/Qfree$ext")
