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
SeparationProblem = CPB.SeparationProblem
Grid = CPB.Grid
empty_grid = CPB.empty_grid

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

βmax = 100.0
N = 2

## Empty
prob = SeparationProblem(N,
                         empty_grid(), empty_grid(),
                         empty_grid(), Grid([[0.0, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep empty" begin
    @test r ≈ 10
    @test norm(af.a, Inf) ≈ 1
end

## Set #1
prob = SeparationProblem(N,
                         Grid([[0.0, 0.0]]), Grid([[0.0, 0.0]]),
                         empty_grid(), Grid([[0.5, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #1.1" begin
    @test r ≈ 0.5/2
    @test norm(af.a, Inf) ≈ 1
end

prob = SeparationProblem(N,
                         Grid([[0.0, 0.0]]), empty_grid(),
                         empty_grid(), Grid([[0.5, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #1.2" begin
    @test r ≈ 0.5
    @test norm(af.a, Inf) ≈ 1
end

## Set #2
prob = SeparationProblem(N,
                         Grid([[0.0, 0.0], [4.0, 0.0]]),
                         Grid([[0.0, 0.0], [4.0, 0.0]]),
                         empty_grid(),
                         Grid([[8.0, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #2.1" begin
    @test r ≈ 4/2
    @test norm(af.a, Inf) ≈ 1
end

prob = SeparationProblem(N,
                         Grid([[4.0, 0.0]]), Grid([[0.0, 0.0]]),
                         empty_grid(), Grid([[8.0, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #2.2" begin
    @test r ≈ 4
    @test norm(af.a, Inf) ≈ 1
end

βmax = 1.0
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #2.3" begin
    @test r ≈ 1
    @test norm(af.a, Inf) ≈ 1
    @test abs(af.β) ≈ βmax
end

## Set #3
βmax = 100.0
prob = SeparationProblem(N,
                         Grid([[1.0, -1.0]]), Grid([[1.0, 1.0]]),
                         empty_grid(), Grid([[2.0, 0.0]]))
af, r = CPB.find_separator(prob, βmax, solver)

@testset "compute sep #3" begin
    @test r ≈ 2/3
    @test norm(af.a, Inf) ≈ 1
end

nothing