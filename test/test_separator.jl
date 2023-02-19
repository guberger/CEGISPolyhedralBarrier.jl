using LinearAlgebra
using JuMP
using HiGHS
using Test
@static if isdefined(Main, :TestLocal) && TestLocal
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

η = 1e-1
βmax = 100.0
N = 2

## Empty
points_inside = Vector{Float64}[]
points_image = Vector{Float64}[]
point_outside = [0.0, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
points_inside = [[0.0, 0.0]]
points_image = [[0.0, 0.0]]
point_outside = [0.5, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #1.1" begin
    @test r ≈ (0.5 - 2*η)/2
    @test norm(af.a, Inf) ≈ 1
end

points_inside = [[0.0, 0.0]]
points_image = Vector{Float64}[]
point_outside = [0.5, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #1.2" begin
    @test r ≈ 0.5/2
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
points_inside = [[0.0, 0.0], [4.0, 0.0]]
points_image = [[0.0, 0.0], [4.0, 0.0]]
point_outside = [8.0, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #2.1" begin
    @test r ≈ (4 - 2*η)/2
    @test norm(af.a, Inf) ≈ 1
end

points_inside = [[4.0, 0.0]]
points_image = [[0.0, 0.0]]
point_outside = [8.0, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #2.2" begin
    @test r ≈ 4/2
    @test norm(af.a, Inf) ≈ 1
end

βmax = 1.0
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #2.3" begin
    @test r ≈ 1/3
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

## Set #3
βmax = 100.0
points_inside = [[1.0, -1.0]]
points_image = [[1.0, 1.0]]
point_outside = [2.0, 0.0]
af, r = CPB.compute_af(
    points_inside, points_image, point_outside, η, βmax, N, solver
)

@testset "compute sep #3" begin
    @test r ≈ (2 - 2*η)/4
    @test norm(af.a, Inf) ≈ 1
end

nothing