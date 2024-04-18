module SpatialOmicsGPs

# Dependencies
using LinearAlgebra
using Statistics
using Optim
using KernelFunctions
using Distances
using Distributions
using DataFrames

# Code files
include("FaST-LMM.jl")
include("SpatialDE.jl")

# FaST-LMM Exports
export fastlmm_fullrank, sigma2_mle_fullrank, minus_log_like_fullrank, delta_mle_fullrank, beta_mle_fullrank_lazy, create_covariate_matrix, project_orth_covar, softplus

# SpatialDE exports
export spatialDE, spatialde_param_inference, spatialde_pvalues, spatialde_varpars_opt, spatialde_varpars_all, spatialde_varpars_null, lengthscales

end
