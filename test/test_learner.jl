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
GenForm = CPB.GenForm
Piece = CPB.Piece
State = CPB.State
BarrierProblem = CPB.BarrierProblem

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

################################################################################
## Set 1

N = 1
pieces = [
    Piece([AffForm([-1.0], 3.0)], 1, [1.0;;], [1.0], 1),
]

prob = BarrierProblem(
    N, pieces,
    [GenForm(1, AffForm([-1.0], 0.0))], # gfs_inv
    [GenForm(1, AffForm([1.0], -4.0))], # gfs_safe
    [State(1, [1.0])], # states_init
    2.0 - eps(2.0), # ϵ
    1e-9 # δ
)

xmax = 1e2
iter_max = 1
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 1: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
    @test length(gen_prob.gfs) == 0
    @test length(gen_prob.states_inside) == 1
    @test length(gen_prob.states_image) == 0
    @test length(gen_prob.links_unknown) == 0
    @test length(gen_prob.links_unknown_new) == 0
    @test length(gen_prob.states_outside) == 0
    @test length(gen_prob.states_outside_new) == 1
end

iter_max = 3
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 1: found" begin
    @test status == CPB.BARRIER_FOUND
    @test length(gen_prob.gfs) == 1
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -1), gen_prob.gfs
    )
    @test length(gen_prob.states_inside) == 1
    @test length(gen_prob.states_image) == 0
    @test length(gen_prob.links_unknown) == 0
    @test length(gen_prob.links_unknown_new) == 0
    @test length(gen_prob.states_outside) == 1
    @test length(gen_prob.states_outside_new) == 0
end

################################################################################
## Set 2

N = 1
pieces = [
    Piece([AffForm([-1.0], 3.0), AffForm([1.0], -3.0)], 1, [1.0;;], [1.5], 1),
    Piece([AffForm([-1.0], 1.0), AffForm([1.0], -1.0)], 1, [0.0;;], [3.5], 1),
]

prob = BarrierProblem(
    N, pieces,
    [GenForm(1, AffForm([-1.0], 0.0))], # gfs_inv
    [GenForm(1, AffForm([1.0], -4.0))], # gfs_safe
    [State(1, [1.0])], # states_init
    1e-4, # ϵ
    1e-9, # δ
)

xmax = 1e2
iter_max = 3
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 2: infeasible" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(gen_prob.states_inside) == 2
    @test length(gen_prob.states_image) == 1
    @test length(gen_prob.states_outside_new) == 1
end

################################################################################
## Set 3

N = 1
pieces = [
    Piece([AffForm([-1.0], 3.0), AffForm([1.0], -3.0)], 1, [0.5;;], [2.5], 1),
    Piece([AffForm([-1.0], 1.0), AffForm([1.0], -1.0)], 2, [0.0;;], [2.5], 1),
    Piece([AffForm([1.0], 0.0)], 2, [0.0;;], [0.9], 2),
]

prob = BarrierProblem(
    N, pieces,
    [GenForm(1, AffForm([-1.0], 0.0))], # gfs_inv
    [GenForm(1, AffForm([1.0], -4.0))], # gfs_safe
    [State(1, [1.0]), State(2, [0.0])], # states_init
    1e2, # ϵ
    1e-3, # δ
)

xmax = 1e2
iter_max = 4
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 3: infeasible" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(gen_prob.states_inside) == 2
    @test length(gen_prob.states_image) == 0
    @test length(gen_prob.states_outside_new) == 1
end

prob = BarrierProblem(
    N, pieces,
    [GenForm(1, AffForm([-1.0], 0.0))], # gfs_inv
    [GenForm(1, AffForm([1.0], -4.0))], # gfs_safe
    [State(1, [1.0]), State(2, [0.0])], # states_init
    1e-2, # ϵ
    1e-3, # δ
)

xmax = 1e2
iter_max = 5
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 3: found 1e-2" begin
    @test status == CPB.BARRIER_FOUND
    @test length(gen_prob.gfs) ≥ 2
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -1), gen_prob.gfs
    )
    @test any(
        gf -> (gf.loc == 2 && gf.af.a ≈ [1] && gf.af.β ≈ -0.95), gen_prob.gfs
    )
end

prob = BarrierProblem(
    N, pieces,
    [GenForm(1, AffForm([-1.0], 0.0))], # gfs_inv
    [GenForm(1, AffForm([1.0], -4.0))], # gfs_safe
    [State(1, [1.0]), State(2, [0.0])], # states_init
    1e-1, # ϵ
    1e-3, # δ
)

xmax = 1e2
iter_max = 5
status, gen_prob = CPB.find_barrier(prob, iter_max, solver, xmax=xmax)

@testset "learn lyapunov 3: found 1e-3" begin
    @test status == CPB.BARRIER_FOUND
    @test length(gen_prob.gfs) ≥ 1
    @test any(
        gf -> (gf.loc == 1 && gf.af.a ≈ [1] && gf.af.β ≈ -2.75), gen_prob.gfs
    )
end

nothing