using SpatialOmicsGPs
using Test

using LinearAlgebra
using Statistics

@testset "SpatialOmicsGPs.jl" begin
    @testset "FaST-LMM" begin
        include("FaST-LMM_tests.jl")
    end
    @testset "SpatialDE" begin
        include("SpatialDE_tests.jl")
    end
end
