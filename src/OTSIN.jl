module OTSIN

using LinearAlgebra

include("optimaltransport.jl")
include("metrics.jl")
include("fitting.jl")
include("utils.jl")
include("perturbation.jl")
include("posterior.jl")

export isprobability, marginals
export optimaltransport
export entropy, relative_entropy, KL, utility
export uniform, report_utility, normalize
export fitM, fitMsgd, r, global_loss
export Regularization, L1, L2, ElasticNet, Exp, Entropic
export perturbation, pert_utility, pert_entropy
export parallel_cross_entropy, parallel_cross_entropy_A, parallel_cross_entropy_B,
        parallel_cross_entropy_free
export laplace_posterior

end # module
