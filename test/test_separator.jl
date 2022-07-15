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

ϵ = 1e-2
βmax = 100.0

## Empty
soft_evids = Point{2}[]
hard_evids = Point{2}[]
point = SVector(0.0, 0.0)
af, r = CPB.compute_af(soft_evids, hard_evids, point, ϵ, βmax, solver)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
soft_evids = [SVector(0.0, 0.0)]
hard_evids = Point{2}[]
point = SVector(0.5, 0.0)
af, r = CPB.compute_af(soft_evids, hard_evids, point, ϵ, βmax, solver)

@testset "compute sep #1" begin
    @test r ≈ 0.5/2 - ϵ/2
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
soft_evids = [SVector(0.0, 0.0)]
hard_evids = [SVector(4.0, 0.0)]
point = SVector(8.0, 0.0)
af, r = CPB.compute_af(soft_evids, hard_evids, point, ϵ, βmax, solver)

@testset "compute sep #2" begin
    @test r ≈ 2 - ϵ
    @test norm(af.a, Inf) ≈ 1
end

## Set #3
βmax = 1.0
af, r = CPB.compute_af(soft_evids, hard_evids, point, ϵ, βmax, solver)

@testset "compute sep #3" begin
    @test r ≈ 1/3 - ϵ
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

nothing