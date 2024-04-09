module SpatialOmicsGPs

# Dependencies
using LinearAlgebra
using Statistics
using Optim
using KernelFunctions
using Distances

# Code files
include("FaST-LMM.jl")


# Exports
export fastlmm_fullrank, sigma2_mle_fullrank, minus_log_like_fullrank, delta_mle_fullrank, project_orth_covar, softplus

end
