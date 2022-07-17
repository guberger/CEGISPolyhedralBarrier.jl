using LinearAlgebra
using StaticArrays
using JuMP
using HiGHS
using Test
@static if isdefined(Main, :TestLocal)
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier
Point = CPB.Point

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

η_inside = 1e-2
η_image = 1e-1
η_outside = 1e-3
βmax = 100.0

## Empty
inside_points = Point{2}[]
image_points = Point{2}[]
outside_point = SVector(0.0, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
inside_points = [SVector(0.0, 0.0)]
image_points = [SVector(0.0, 0.0)]
outside_point = SVector(0.5, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #1.1" begin
    @test r ≈ (0.5 - η_image - η_outside)/2
    @test norm(af.a, Inf) ≈ 1
end

inside_points = [SVector(0.0, 0.0)]
image_points = Point{2}[]
outside_point = SVector(0.5, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #1.2" begin
    @test r ≈ (0.5 - η_inside - η_outside)/2
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
inside_points = [SVector(0.0, 0.0), SVector(4.0, 0.0)]
image_points = [SVector(0.0, 0.0), SVector(4.0, 0.0)]
outside_point = SVector(8.0, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #2.1" begin
    @test r ≈ (4 - η_image - η_outside)/2
    @test norm(af.a, Inf) ≈ 1
end

inside_points = [SVector(4.0, 0.0)]
image_points = [SVector(0.0, 0.0)]
outside_point = SVector(8.0, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #2.2" begin
    @test r ≈ (4 - η_inside - η_outside)/2
    @test norm(af.a, Inf) ≈ 1
end

## Set #3
βmax = 1.0
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #3" begin
    @test r ≈ (1 - 2*η_inside - η_outside)/3
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

## Set #4
βmax = 100.0
inside_points = [SVector(1.0, -1.0)]
image_points = [SVector(1.0, 1.0)]
outside_point = SVector(2.0, 0.0)
af, r = CPB.compute_af(
    inside_points, image_points, outside_point,
    η_inside, η_image, η_outside, βmax, solver
)

@testset "compute sep #4" begin
    @test r ≈ (2 - η_inside - η_image - 2*η_outside)/4
    @test norm(af.a, Inf) ≈ 1
end

nothing