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
State = CPB.State
Piece = CPB.Piece
VerifierProblem = CPB.VerifierProblem
CexKey = CPB.CexKey
CexVal = CPB.CexVal

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

xmax = 1e3

################################################################################
## Set #1

N = 1
pieces = [Piece([AffForm([-1.0], 0.0)], 1, [0.5;;], [1.0], 1)]

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_bf
    GenForm[], # gfs_safe
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: empty" begin
    @test isempty(prob.keys_todo)
    @test isempty(prob.cexs)
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_bf
    [GenForm(1, AffForm([1.0], -0.25))], # gfs_safe
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: r neg" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [0.375]
    @test prob.cexs[CexKey(1, 1)].r ≈ -0.125
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_bf
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_safe
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: r pos" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [0.75]
    @test prob.cexs[CexKey(1, 1)].r ≈ 0.25
end

empty!(prob.keys_todo)
prob.cexs[CexKey(1, 1)].state.x[1] = -5.0
CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: add infeasible" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [0.75]
    @test prob.cexs[CexKey(1, 1)].r ≈ 0.25
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_bf
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_safe
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_bf, eachindex(prob.gfs_bf))
CPB.update_cexs_cont!(prob, xmax, solver)

@testset "verify cont: r neg" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test norm(prob.cexs[CexKey(1, 1)].state.x) < 1e-6
    @test prob.cexs[CexKey(1, 1)].r ≈ -0.5
end

empty!(prob.keys_todo)
prob.cexs[CexKey(1, 1)].state.x[1] = -5.0
CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.update_cexs_cont!(prob, xmax, solver)

@testset "verify cont: add infeasible" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test norm(prob.cexs[CexKey(1, 1)].state.x) < 1e-6
    @test prob.cexs[CexKey(1, 1)].r ≈ -0.5
end

################################################################################
## Set 2

N = 2
pieces = [
    Piece([AffForm([-1.0, 0.0], -1.0)], 1, [0.0 1.0; -1.0 0.0], [0.0, 0.0], 2)
]

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0, 0.0], -1.0))], # gfs_bf
    [GenForm(1, AffForm([0.0, -1.0], -1.0))], # gfs_safe
    [GenForm(1, AffForm([0.0, 1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: empty" begin
    @test isempty(prob.keys_todo)
    @test isempty(prob.cexs)
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0, 0.0], -1.0))], # gfs_bf
    [GenForm(1, AffForm([0.0, -1.0], -1.0))], # gfs_safe
    [GenForm(1, AffForm([0.0, 1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_bf, eachindex(prob.gfs_bf))
CPB.update_cexs_cont!(prob, xmax, solver)

@testset "verify cont: empty" begin
    @test isempty(prob.keys_todo)
    @test isempty(prob.cexs)
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(2, AffForm([1.0, 0.0], -1.0))], # gfs_bf
    [GenForm(1, AffForm([0.0, -1.0], -1.0))], # gfs_safe
    [GenForm(1, AffForm([0.0, 1.0], -1.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_bf, eachindex(prob.gfs_bf))
CPB.update_cexs_cont!(prob, xmax, solver)

@testset "verify cont: r zero" begin
    @test length(prob.keys_todo) == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x[2] ≈ 1
    @test abs(prob.cexs[CexKey(1, 1)].r) < 1e-6
end

################################################################################
## Set 3

N = 2
pieces = [
    Piece([AffForm([-1.0, 0.0], -1.0)], 1, [0.0 1.0; -1.0 0.0], [0.0, +0.0], 2),
    Piece([AffForm([-1.0, 0.0], -1.0)], 1, [1.0 0.0; 0.0 1.0], [0.0, -1.0], 2)
]

prob = VerifierProblem(
    N, pieces,
    GenForm[], # gfs_bf
    [
        GenForm(1, AffForm([0.0, -1.0], -1.0)),
        GenForm(2, AffForm([0.0, -1.0], -1.0))
    ], # gfs_safe
    [GenForm(1, AffForm([1.0, 0.0], -2.0))], # gfs_inv
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:4]..., # all afs_...
)

CPB.add_infeasible_keys!(prob, prob.gfs_bf)
CPB.add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
CPB.update_cexs_safe!(prob, xmax, solver)

@testset "verify safe: r pos" begin
    @test length(prob.keys_todo) == 2
    @test length(prob.cexs) == 2
    @test prob.cexs[CexKey(1, 2)].state.loc == 1
    @test prob.cexs[CexKey(1, 2)].state.x[1] ≈ 1
    @test prob.cexs[CexKey(1, 2)].r ≈ xmax + 1
    @test prob.cexs[CexKey(2, 2)].state.loc == 1
    @test prob.cexs[CexKey(2, 2)].state.x[2] ≈ 0
    @test prob.cexs[CexKey(2, 2)].r ≈ 1
end

nothing