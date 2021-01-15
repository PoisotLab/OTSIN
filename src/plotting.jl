#=
Created on Tuesday 14 July 2020
Last update: -

@author: Michiel Stock
michielfmstock@gmail.com

Plotting utility for OTSIN
=#

using Plots, RecipesBase, LaTeXStrings


@userplot MarginalsHeatmap

@recipe function f(h::MarginalsHeatmap; cmap=cgrad([:white, :blue]), cola="#2a9d8f", colb="#f4a261", showaxis=false)
    if length(h.args) != 1 || !(typeof(h.args[1]) <: AbstractMatrix)
        error("Marginals Heatmap should be given a matrix.  Got: $(typeof(h.args))")
    end
    P, = h.args
    a, b = marginals(P)

    # set up the subplots
    legend := false
    link := :both
    framestyle := [:none :axes :none]
    grid := false
    layout := @layout [tophist           _
                       hist2d{0.9w,0.9h} righthist]

    # main histogram2d
    @series begin
        seriescolor := cmap
        seriestype := :heatmap
        subplot := 2
        showaxis := showaxis
        P
    end

    # these are common to both marginal histograms
    fillalpha := 0.8
    linealpha := 0.3
    seriestype := :bar

    # upper histogram
    @series begin
        fillcolor := cola
        subplot := 3
        orientation := :h
        title := L"$\mathbf{a}$"
        a
    end

    # right histogram
    @series begin
        fillcolor := colb
        subplot := 1
        title := L"$\mathbf{b}$"
        b
    end
end
