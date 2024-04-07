using SpatialOmicsGPs
using Test

using LinearAlgebra
using Statistics

@testset "SpatialTranscriptomicsGPs.jl" begin
    @testset "FaSTLMM" begin
        @testset "FaSTLMM" begin
            include("FaST-LMM_tests.jl")
        end
    end
end
