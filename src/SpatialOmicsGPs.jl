module SpatialOmicsGPs

# Dependencies
using LinearAlgebra
using Statistics
#using Optim
#using KernelFunctions
#using Distances
using Distributions
using DataFrames
using FaSTLMMlight

# Code files
include("SpatialDE.jl")

# SpatialDE exports
export spatialDE, spatialde_param_inference, spatialde_pvalues, spatialde_varpars_opt, spatialde_varpars_all, spatialde_varpars_null, lengthscales

end
