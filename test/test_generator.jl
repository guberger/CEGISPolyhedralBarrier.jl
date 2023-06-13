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
AffForm = CPB.AffForm
PolyFunc = CPB.PolyFunc
empty_pf = CPB.empty_pf
Grid = CPB.Grid
empty_grid = CPB.empty_grid
Link = CPB.Link
Graph = CPB.Graph
empty_graph = CPB.empty_graph
GeneratorProblem = CPB.GeneratorProblem

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

N = 1
ϵ = 1.5
βmax = 100.0

################################################################################
## Generator: feasible

prob = GeneratorProblem(
    N,
    [PolyFunc([AffForm([0.0], 0.0)]) for loc = 1:2], # mpf
    [empty_pf() for loc = 1:2], # mpf_new
    [true for loc = 1:2], # mreset
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
    @test all(prob.mreset)
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

################################################################################
## Generator: infeasible

prob = GeneratorProblem(
    N,
    [PolyFunc([AffForm([0.0], 0.0)]) for loc = 1:2], # mpf
    [empty_pf() for loc = 1:2], # mpf_new
    [true for loc = 1:2], # mreset
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
    @test all(prob.mreset)
    @test length(prob.mgrid_inside[1].points) == 2
    @test length(prob.mgrid_inside[2].points) ≥ 1
    @test length(prob.mgrid_image[1].points) == 1
    @test length(prob.mgrid_image[2].points) ≥ 1
end

################################################################################
## Generator: partially reset

prob = GeneratorProblem(
    N,
    [PolyFunc([AffForm([0.0], 0.0)]) for loc = 1:2], # mpf
    [PolyFunc([AffForm([0.0], 0.0)]) for loc = 1:2], # mpf_new
    [true for loc = 1:2], # mreset
    [Grid([[0.0]]), empty_grid()], # mgrid_inside
    [Grid([[0.0]]) for loc = 1:2], # mgrid_image
    Graph([
        Link(1, [0.0], 2, [9.0]),
        Link(2, [9.0], 1, [9.0])
    ]), # graph_unknown
    Graph([
        Link(1, [5.0], 2, [2.0]),
        Link(2, [2.0], 2, [2.0])
    ]), # graph_unknown_new
    Graph([
        Link(1, [0.0], 0, [NaN]),
        Link(2, [-9.0], 0, [NaN])
    ]), # graph_outside
    Graph([
        Link(1, [-5.0], 0, [NaN])
    ]), # graph_outside_new
    empty_graph(), # graph_temp
    ϵ
)

flag = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: infeasible" begin
    @test flag
    @test !prob.mreset[1]
    @test prob.mreset[2]
    @test length(prob.mpf[1].afs) == 3
    @test length(prob.mpf[2].afs) == 2
    @test length(prob.mpf_new[1].afs) == 2
    @test length(prob.mgrid_inside[1].points) == 1
    @test length(prob.mgrid_inside[2].points) == 1
    @test length(prob.mgrid_image[1].points) == 1
    @test length(prob.mgrid_image[2].points) == 2
    @test length(prob.graph_outside.links) == 3
    @test length(prob.graph_unknown.links) == 3
    for pf in (prob.mpf[1], prob.mpf_new[1])
        @test any(
            af -> (af.a ≈ [1] && af.β ≈ -2.5), pf.afs
        )
        @test any(
            af -> (af.a ≈ [-1] && af.β ≈ -2.5), pf.afs
        )
    end
    @test any(
        af -> (af.a ≈ [1] && af.β ≈ -5.5), prob.mpf[2].afs
    )
    @test any(
        af -> (af.a ≈ [-1] && af.β ≈ -4.5), prob.mpf[2].afs
    )
end

nothing