using StatsBase

function perturbation(f, M, a, b; λ::Real=1.0, maxitter::Int=10, ϵ::Real=1e-3)
        feval = (a, b) -> f(optimaltransport(M, a, b, λ=λ, maxitter=maxitter, ϵ=ϵ))
        ∇af, ∇bf = Zygote.gradient(feval, a, b)
        return ∇af .- mean(∇af), ∇bf .- mean(∇bf)
end


function pert_utility(M, a, b; λ::Real=1.0, maxitter::Int=10, ϵ::Real=1e-3)
   return perturbation(Q -> utility(Q, M), M, a, b, λ=λ, maxitter=maxitter, ϵ=ϵ)
end

function pert_entropy(M, a, b; λ::Real=1.0, maxitter::Int=10, ϵ::Real=1e-3)
   return perturbation(Q -> entropy(Q), M, a, b, λ=λ, maxitter=maxitter, ϵ=ϵ)
end
