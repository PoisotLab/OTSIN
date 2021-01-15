### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9e3f6282-3575-11eb-35e8-69e4cf57b6f4
using LinearAlgebra, Plots, PlutoUI

# ╔═╡ 16df5992-3574-11eb-3682-0981f8882895
md"""
# Generalized optimal transport

Traditional optimal transport find the optimal transportation map given fixed marginals. For many applications, this requirement is too stringent, for example, because the marginals are not balanced. So, one can define **generalized optimal transport** as:

$$\max_Q\, \langle M, Q \rangle + \epsilon H(Q) + \alpha D_\text{KL}(\mathbf{a}\mid Q \mathbb{1}) + \beta D_\text{KL}(\mathbf{b}\mid Q^\intercal \mathbb{1})\,.$$

One notes that the **Kullback-Leibler** is now used as a soft-constraints to drive the marginals of $Q$ to match $\mathbf{a}$ and $\mathbf{b}$, these themselfs don't have to be normalized.

One can solve this problem using the **generalized Sinkhorn algorithm**.

"""

# ╔═╡ 244cb164-3575-11eb-3b9d-837acaa0e9e9
function gen_sinkhorn(M, a, b; ϵ=1.0, α=1.0, β=1.0)
	K = exp.(M ./ ϵ)
	u, v = ones(length(a)), ones(length(b))
	for i in 1:10
		u .= (a ./ K * v) .^ (α / (α + ϵ))
		v .= (b ./ K' * u) .^ (β / (β + ϵ))
	end
	K .*= u
	K .*= v'
	K ./= sum(K)
	return K
end

# ╔═╡ 62aa2bf2-3576-11eb-0cc6-e31f7ad771d0
n, m = 20, 10

# ╔═╡ 3486c424-3576-11eb-35a3-0bcc30b65445
M = randn(n, m)

# ╔═╡ ac95ae3a-357b-11eb-1e16-b56d76851314
normalize(x) = x ./ sum(x)

# ╔═╡ 5ecc17ca-3576-11eb-00ce-2fe7c5d8f55b
a = rand(n) #|> normalize

# ╔═╡ 72fc2c12-3576-11eb-23d1-c1d14da6de89
b = rand(m) #|> normalize

# ╔═╡ 3d0844ac-3577-11eb-3086-112eb6425f79
plot(bar(a), bar(b), layout=(2,1))

# ╔═╡ d9ffdbac-3576-11eb-24f1-55d257632108
@bind α Slider(0:0.01:10, default=1)

# ╔═╡ f22bc95c-3576-11eb-2a14-9d92a4497f68
@bind β Slider(0:0.01:10, default=1)

# ╔═╡ cd231944-357b-11eb-24b6-dbbb08cb40d3
@bind ϵ Slider(0:0.01:10, default=1)

# ╔═╡ 7702f944-3576-11eb-08e4-0f1f61b1d3bd
Q = gen_sinkhorn(M, a, b; ϵ, α, β)

# ╔═╡ c80cf7c2-3576-11eb-0e81-05121d909271
sum(Q)

# ╔═╡ d0a700bc-3576-11eb-2d6c-2354903c6a7d
heatmap(Q)

# ╔═╡ 84a3ccea-357b-11eb-0924-b9f6fe8c481f
ã = sum(Q, dims=2)[:]

# ╔═╡ 8fd9f8dc-357b-11eb-0fa8-812f0ea106e5
b̃ = sum(Q, dims=1)[:]

# ╔═╡ 7c05eabe-357b-11eb-0f67-4fa981d4d725
plot(bar(ã), bar(b̃), layout=(2,1))

# ╔═╡ Cell order:
# ╠═16df5992-3574-11eb-3682-0981f8882895
# ╠═9e3f6282-3575-11eb-35e8-69e4cf57b6f4
# ╠═244cb164-3575-11eb-3b9d-837acaa0e9e9
# ╠═62aa2bf2-3576-11eb-0cc6-e31f7ad771d0
# ╠═3486c424-3576-11eb-35a3-0bcc30b65445
# ╠═ac95ae3a-357b-11eb-1e16-b56d76851314
# ╠═5ecc17ca-3576-11eb-00ce-2fe7c5d8f55b
# ╠═72fc2c12-3576-11eb-23d1-c1d14da6de89
# ╠═3d0844ac-3577-11eb-3086-112eb6425f79
# ╠═d9ffdbac-3576-11eb-24f1-55d257632108
# ╠═f22bc95c-3576-11eb-2a14-9d92a4497f68
# ╠═cd231944-357b-11eb-24b6-dbbb08cb40d3
# ╠═7702f944-3576-11eb-08e4-0f1f61b1d3bd
# ╠═c80cf7c2-3576-11eb-0e81-05121d909271
# ╠═d0a700bc-3576-11eb-2d6c-2354903c6a7d
# ╠═84a3ccea-357b-11eb-0924-b9f6fe8c481f
# ╠═8fd9f8dc-357b-11eb-0fa8-812f0ea106e5
# ╠═7c05eabe-357b-11eb-0f67-4fa981d4d725
