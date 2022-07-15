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

η = 1e-2
βmax = 100.0

## Empty
neg_points = Point{2}[]
point = SVector(0.0, 0.0)
af, r = CPB.compute_af(neg_points, point, η, βmax, solver)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
neg_points = [SVector(0.0, 0.0)]
point = SVector(0.5, 0.0)
af, r = CPB.compute_af(neg_points, point, η, βmax, solver)

@testset "compute sep #1" begin
    @test r ≈ 0.5/2 - η
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
neg_points = [SVector(0.0, 0.0), SVector(4.0, 0.0)]
point = SVector(8.0, 0.0)
af, r = CPB.compute_af(neg_points, point, η, βmax, solver)

@testset "compute sep #2" begin
    @test r ≈ 2 - η
    @test norm(af.a, Inf) ≈ 1
end

## Set #3
βmax = 1.0
af, r = CPB.compute_af(neg_points, point, η, βmax, solver)

@testset "compute sep #3" begin
    @test r ≈ 1/3 - η
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

nothing