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
GeneratorProblem = CPB.GeneratorProblem
PolyFunc = CPB.PolyFunc
empty_pf = CPB.empty_pf
Grid = CPB.Grid
empty_grid = CPB.empty_grid
Link = CPB.Link
Graph = CPB.Graph
empty_graph = CPB.empty_graph

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

N = 1
ϵ = 1.5
βmax = 100.0

prob = GeneratorProblem(
    N,
    [empty_pf() for loc = 1:2], # mpf
    [Grid([[0.0]]), empty_grid()], # mgrid_inside
    [empty_grid() for loc = 1:2], # mgrid_image
    Graph([
        Link(2, [4.0], 1, [4.0]),
        Link(1, [9.0], 2, [9.0]),
        Link(2, [9.0], 1, [9.0])
    ]), # graph_unknown
    Graph([
        Link(1, [1.0], 2, [2.0])
    ]), # graph_unknown_new
    Graph([
        Link(1, [-4.0], 0, [NaN]),
    ]), # graph_outside
    empty_graph(), # graph_outside_new
    empty_graph(), # graph_temp
    ϵ
)

flag = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: feasible" begin
    @test flag
    @test length(prob.mpf[1].afs) == 2
    @test length(prob.mpf[2].afs) == 1
    @test length(prob.mgrid_inside[1].points) == 2
    @test length(prob.mgrid_inside[2].points) == 1
    @test length(prob.mgrid_image[1].points) == 1
    @test length(prob.mgrid_image[2].points) == 1
    @test length(prob.graph_outside.links) == 1
    @test length(prob.graph_unknown.links) == 2
    @test any(
        af -> (af.a ≈ [1] && af.β ≈ -6.5), prob.mpf[1].afs
    )
    @test any(
        af -> (af.a ≈ [-1] && abs(af.β) < 1e-6), prob.mpf[1].afs
    )
    @test any(
        af -> (af.a ≈ [1] && af.β ≈ -5.5), prob.mpf[2].afs
    )
end

prob = GeneratorProblem(
    N,
    [empty_pf() for loc = 1:2], # mpf
    [Grid([[0.0]]), empty_grid()], # mgrid_inside
    [empty_grid() for loc = 1:2], # mgrid_image
    Graph([
        Link(2, [4.0], 1, [4.0]),
        Link(1, [9.0], 2, [9.0]),
        Link(2, [9.0], 1, [9.0])
    ]), # graph_unknown
    Graph([
        Link(1, [1.0], 2, [2.0])
    ]), # graph_unknown_new
    Graph([
        Link(1, [-4.0], 0, [NaN]),
    ]), # graph_outside
    Graph([
        Link(1, [5.0], 0, [NaN]),
    ]), # graph_outside_new
    empty_graph(), # graph_temp
    ϵ
)

flag = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: infeasible" begin
    @test !flag
    @test length(prob.mgrid_inside[1].points) == 2
    @test length(prob.mgrid_inside[2].points) ≥ 1
    @test length(prob.mgrid_image[1].points) == 1
    @test length(prob.mgrid_image[2].points) ≥ 1
end

nothing