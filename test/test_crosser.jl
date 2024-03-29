using LinearAlgebra
using JuMP
using HiGHS
using Test
@static if isdefined(Main, :TestLocal)
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
CrossingProblem = CPB.CrossingProblem

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

xmax = 1e3

# Set #1
N = 1
A = [0.5;;]
b = [1.0]

prob = CrossingProblem(
    N, A, b,
    [AffForm([-1.0], 0.0), AffForm([1.0], 1.0), AffForm([-0.5], 0.25)],
    AffForm[]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser safe: infeasible" begin
    @test r == -Inf
    @test all(isnan, x)
    @test !flag
end

prob = CrossingProblem(
    N, A, b,
    [AffForm([-1.0], 0.0), AffForm([1.0], -1.0), AffForm([-0.5], 0.25)],
    AffForm[]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser safe: empty outside" begin
    @test r ≈ 10*xmax
    @test flag
end

prob = CrossingProblem(
    N, A, b,
    [AffForm([1.0], -1.0), AffForm([-0.5], 0.25)],
    [AffForm([1.0], -1.0)]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser safe: positive radius" begin
    @test r ≈ 0.5
    @test x ≈ [1]
    @test flag
end

prob = CrossingProblem(
    N, A, b,
    [AffForm([1.0], -1.0), AffForm([-0.5], 0.25)],
    [AffForm([-0.5], 0.25)]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser BF: safe" begin
    @test r ≈ -0.75
    @test x ≈ [0.5]
    @test flag
end

# Set #2
N = 2
A = [0.0 1.0; -1.0 0.0]
b = [0.0, 0.0]

prob = CrossingProblem(
    N, A, b,
    [
        AffForm([-1.0, 0.0], -1.0), AffForm([0.0, 1.0], -1.0),
        AffForm([0.0, -1.0], -1.0), AffForm([1.0, 0.0], -1.0)
    ],
    AffForm[]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser safe: empty outside" begin
    @test r ≈ 10*xmax
    @test flag
end

prob = CrossingProblem(
    N, A, b,
    [
        AffForm([-1.0, 0.0], -1.0),
        AffForm([0.0, 1.0], -1.0),
        AffForm([0.0, -1.0], -1.0)
    ],
    [AffForm([1.0, 0.0], -1.0)]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser BF: just safe" begin
    @test abs(r) < 1e-6
    @test x[2] ≈ 1
    @test flag
end

prob = CrossingProblem(
    N, A, b,
    [
        AffForm([-1.0, 0.0], -1.0),
        AffForm([0.0, 1.0], -1.0),
        AffForm([0.0, -1.0], -1.0)
    ],
    [AffForm([0.0, -1.0], -1.0)]
)

x, r, flag = CPB.find_crosser(prob, xmax, solver)

@testset "crosser safe: totally unsafe" begin
    @test r ≈ xmax - 1
    @test x[1] ≈ xmax
    @test flag
end

nothing