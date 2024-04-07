using SpatialOmicsGPs
using Test

using LinearAlgebra
using Statistics

@testset "SpatialOmicsGPs.jl" begin
    @testset "FaSTLMM" begin
        @testset "FaST-LMM" begin
            include("FaST-LMM_tests.jl")
        end
    end
end
