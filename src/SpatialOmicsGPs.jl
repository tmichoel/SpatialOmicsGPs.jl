module SpatialTranscriptomicsGPs

# Dependencies
using LinearAlgebra
using Statistics
using Optim

# Code files
include("FaSTLMM.jl")


# Exports
export beta_mle_fullrank, sigma2_mle_fullrank, minus_log_like_fullrank, delta_mle_fullrank

end
