using SpatialOmicsGPs
using Test

using LinearAlgebra
using Statistics
using FaSTLMMlight

@testset "SpatialOmicsGPs.jl" begin
    @testset "SpatialDE" begin
        include("SpatialDE_tests.jl")
    end
end
