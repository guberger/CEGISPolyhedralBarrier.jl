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
GenForm = CPB.GenForm
State = CPB.State
Link = CPB.Link
GeneratorProblem = CPB.GeneratorProblem
empty_xs() = Vector{Float64}[]

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
    [GenForm(0, AffForm([0.0], 0.0))], # gfs
    [GenForm(1, AffForm([1.0], -10.0))], # gfs_safe
    [0], # indices_new
    [State(1, [1.0])], # states_inside
    State[], # states_image
    [
        Link(State(2, [4.0]), State(1, [4.0])),
        Link(State(1, [9.0]), State(2, [9.0])),
        Link(State(2, [9.0]), State(1, [9.0]))
    ], # links_unknown
    [
        Link(State(1, [1.0]), State(2, [2.0]))
    ], # links_unknown_new
    [State(1, [-3.0])], # states_outside
    State[], # states_outside_new
    [empty_xs() for i = 1:4]..., # all xs_...
    ϵ
)

isreset, issuccess = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: feasible" begin
    @test isreset
    @test issuccess
    @test length(prob.gfs) == 4
    @test prob.indices_new == collect(eachindex(prob.gfs))
    @test length(prob.states_inside) == 3
    @test length(prob.states_image) == 2
    @test length(prob.links_unknown) == 2
    @test length(prob.states_outside) == 1
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -7), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -6.5), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [-1] && gf.af.β ≈ 0.5), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 2 && gf.af.a ≈ [1] && gf.af.β ≈ -5.5), prob.gfs
    )
end

################################################################################
## Generator: infeasible (point outside)

prob = GeneratorProblem(
    N,
    [GenForm(0, AffForm([0.0], 0.0))], # gfs
    [GenForm(1, AffForm([0.0], -10.0))], # gfs_safe
    [0], # indices_new
    [State(1, [1.0])], # states_inside
    State[], # states_image
    [
        Link(State(2, [4.0]), State(1, [4.0])),
        Link(State(1, [9.0]), State(2, [9.0])),
        Link(State(2, [9.0]), State(1, [9.0]))
    ], # links_unknown
    [
        Link(State(1, [1.0]), State(2, [2.0]))
    ], # links_unknown_new
    [State(1, [-3.0])], # states_outside
    [State(1, [5.0])], # states_outside_new
    [empty_xs() for i = 1:4]..., # all xs_...
    ϵ
)

isreset, issuccess = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: infeasible (point outside)" begin
    @test isreset
    @test !issuccess
    @test prob.indices_new == collect(eachindex(prob.gfs))
    @test length(prob.states_inside) == 3
    @test length(prob.states_image) == 2
end

################################################################################
## Generator: infeasible (safe)

prob = GeneratorProblem(
    N,
    [GenForm(0, AffForm([0.0], 0.0))], # gfs
    [GenForm(1, AffForm([1.0], -5.0))], # gfs_safe
    [0], # indices_new
    [State(1, [1.0])], # states_inside
    State[], # states_image
    [
        Link(State(2, [4.0]), State(1, [4.0])),
        Link(State(1, [9.0]), State(2, [9.0])),
        Link(State(2, [9.0]), State(1, [9.0]))
    ], # links_unknown
    [
        Link(State(1, [1.0]), State(2, [2.0]))
    ], # links_unknown_new
    [State(1, [-3.0])], # states_outside
    State[], # states_outside_new
    [empty_xs() for i = 1:4]..., # all xs_...
    ϵ
)

isreset, issuccess = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: infeasible (safe)" begin
    @test isreset
    @test !issuccess
    @test prob.indices_new == collect(eachindex(prob.gfs))
    @test length(prob.states_inside) == 3
    @test length(prob.states_image) == 2
end

################################################################################
## Generator: no reset

prob = GeneratorProblem(
    N,
    [GenForm(0, AffForm([0.0], 0.0))], # gfs
    [GenForm(1, AffForm([1.0], -10.0))], # gfs_safe
    [0], # indices_new
    State[], # states_inside
    [State(1, [1.0])], # states_image
    [
        Link(State(2, [4.0]), State(1, [4.0])),
        Link(State(1, [9.0]), State(2, [9.0])),
        Link(State(2, [9.0]), State(1, [9.0]))
    ], # links_unknown
    [
        Link(State(1, [-5.0]), State(2, [2.0]))
    ], # links_unknown_new
    [State(1, [1.0])], # states_outside
    [State(1, [5.0])], # states_outside_new
    [empty_xs() for i = 1:4]..., # all xs_...
    ϵ
)

isreset, issuccess = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: no reset" begin
    @test !isreset
    @test issuccess
    @test length(prob.gfs) == 3
    @test prob.indices_new == collect(2:length(prob.gfs))
    @test length(prob.states_inside) == 0
    @test length(prob.states_image) == 1
    @test length(prob.links_unknown) == 4
    @test length(prob.states_outside) == 2
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [-1] && gf.af.β ≈ -2), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -3), prob.gfs
    )
end

prob = GeneratorProblem(
    N,
    State[], # gfs
    [GenForm(1, AffForm([1.0], -10.0))], # gfs_safe
    [0], # indices_new
    State[], # states_inside
    [State(1, [1.0])], # states_image
    [
        Link(State(2, [4.0]), State(1, [4.0])),
        Link(State(1, [9.0]), State(2, [9.0])),
        Link(State(2, [9.0]), State(1, [9.0]))
    ], # links_unknown
    [
        Link(State(1, [-5.0]), State(2, [2.0]))
    ], # links_unknown_new
    [State(1, [1.0])], # states_outside
    [State(1, [5.0])], # states_outside_new
    [empty_xs() for i = 1:4]..., # all xs_...
    ϵ
)

isreset, issuccess = CPB.update_generator!(prob, βmax, solver)

@testset "Generator: starts empty" begin
    @test isreset
    @test issuccess
    @test length(prob.gfs) == 3
    @test prob.indices_new == collect(eachindex(prob.gfs))
    @test length(prob.states_inside) == 0
    @test length(prob.states_image) == 1
    @test length(prob.links_unknown) == 4
    @test length(prob.states_outside) == 2
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -5.5), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [-1] && gf.af.β ≈ -2), prob.gfs
    )
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -3), prob.gfs
    )
end

nothing