using SpatialTranscriptomicsGPs
using Test

using LinearAlgebra
using Statistics

@testset "SpatialTranscriptomicsGPs.jl" begin
    @testset "FaSTLMM" begin
        @testset "FaSTLMM" begin
            include("FaSTLMM_tests.jl")
        end
    end
end
