module OTSIN

using LinearAlgebra

include("optimaltransport.jl")
include("metrics.jl")
include("fitting.jl")
include("utils.jl")
include("perturbation.jl")

export isprobability, marginals
export optimaltransport
export entropy, relative_entropy, KL, utility
export uniform, report_utility, normalize
export fitM, fitMsgd, r
export Regularization, L1, L2, ElasticNet, Exp, Entropic
export perturbation, pert_utility, pert_entropy

end # module
