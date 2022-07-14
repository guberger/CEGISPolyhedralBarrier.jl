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
Separator = CPB.Separator

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

ϵ = 1e-2
βmax = 100.0

## Empty
sep = Separator{2}(ϵ, βmax)
point = SVector(0.0, 0.0)
af, r = CPB.compute_af(sep, point, solver)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
sep = Separator{2}(ϵ, βmax)
CPB.add_soft_evid!(sep, SVector(0.0, 0.0))
point = SVector(0.5, 0.0)
af, r = CPB.compute_af(sep, point, solver)

@testset "compute sep #1" begin
    @test r ≈ 0.5/2 - ϵ/2
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
sep = Separator{2}(ϵ, βmax)
CPB.add_soft_evid!(sep, SVector(0.0, 0.0))
CPB.add_hard_evid!(sep, SVector(4.0, 0.0))
point = SVector(8.0, 0.0)
af, r = CPB.compute_af(sep, point, solver)

@testset "compute sep #2" begin
    @test r ≈ 2 - ϵ
    @test norm(af.a, Inf) ≈ 1
end

## Set #3
βmax = 1.0
sep = Separator{2}(ϵ, βmax)
CPB.add_soft_evid!(sep, SVector(0.0, 0.0))
CPB.add_hard_evid!(sep, SVector(4.0, 0.0))
point = SVector(8.0, 0.0)
af, r = CPB.compute_af(sep, point, solver)

@testset "compute sep #3" begin
    @test r ≈ 1/3 - ϵ
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

nothing