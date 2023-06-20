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
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_bf
    GenForm[], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier: empty" begin
    @test nkey == 0
    @test isempty(prob.cexs)
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    [
        GenForm(1, AffForm([1.0], -0.25)),
        GenForm(1, AffForm([-0.5], 0.25))
    ], # gfs_bf
    [GenForm(1, AffForm([1.0], -0.25))], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier: infeasible" begin
    @test nkey == 1
    @test isempty(prob.cexs)
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    [
        GenForm(1, AffForm([1.0], -1.0)),
        GenForm(1, AffForm([-0.5], 0.25))
    ], # gfs_bf
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier safe: r pos" begin
    @test nkey == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [1]
    @test prob.cexs[CexKey(1, 1)].r ≈ 0.5
end

empty!(prob.keys_todo)
prob.cexs[CexKey(1, 1)].state.x[1] = 0.5
CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier safe: add infeasible nothing" begin
    @test nkey == 0
    @test length(prob.cexs) == 1
end

empty!(prob.keys_todo)
prob.cexs[CexKey(1, 1)].state.x[1] = 0.25
CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier safe: add infeasible" begin
    @test nkey == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [1]
    @test prob.cexs[CexKey(1, 1)].r ≈ 0.5
end

#-------------------------------------------------------------------------------
prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0], -1.0))], # gfs_inv
    [
        GenForm(1, AffForm([1.0], -1.0)),
        GenForm(1, AffForm([-0.5], 0.25))
    ], # gfs_bf
    [GenForm(1, AffForm([-0.5], 0.25))], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier: r neg" begin
    @test nkey == 1
    @test length(prob.cexs) == 1
    @test prob.cexs[CexKey(1, 1)].state.loc == 1
    @test prob.cexs[CexKey(1, 1)].state.x ≈ [0.5]
    @test prob.cexs[CexKey(1, 1)].r ≈ -0.75
end

################################################################################
## Set 2

N = 2
pieces = [
    Piece([AffForm([-1.0, 0.0], -1.0)], 1, [0.0 1.0; -1.0 0.0], [0.0, 0.0], 2)
]

prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([0.0, 1.0], -1.0))], # gfs_inv
    [
        GenForm(1, AffForm([0.0, -1.0], -1.0)),
        GenForm(2, AffForm([1.0, 0.0], -1.0))
    ], # gfs_bf
    [GenForm(2, AffForm([1.0, 0.0], -1.0))], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier: r zero" begin
    @test nkey == 1
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
    Piece([AffForm([-1.0, 0.0], -1.0)], 1, [1.0 0.0; 0.0 1.0], [0.0, -2.0], 2)
]

prob = VerifierProblem(
    N, pieces,
    [GenForm(1, AffForm([1.0, 0.0], -2.0))], # gfs_inv
    [
        GenForm(1, AffForm([0.0, -1.0], -1.0)),
        GenForm(2, AffForm([0.0, -1.0], -1.0))
    ], # gfs_bf
    [
        GenForm(1, AffForm([0.0, -1.0], -1.0)),
        GenForm(2, AffForm([0.0, -1.0], -1.0))
    ], # gfs_out
    Dict{CexKey,CexVal}(),
    CexKey[],
    [AffForm[] for i = 1:2]... # all afs_...
)

CPB.add_keys_bf_infeasible!(prob, eachindex(prob.gfs_bf))
CPB.add_keys_out_new!(prob, eachindex(prob.gfs_out))
nkey = length(prob.keys_todo)
CPB.update_cexs!(prob, xmax, solver)

@testset "Verifier: r pos" begin
    @test nkey == 2
    @test length(prob.cexs) == 2
    @test prob.cexs[CexKey(1, 2)].state.loc == 1
    @test prob.cexs[CexKey(1, 2)].state.x[1] ≈ 2
    @test prob.cexs[CexKey(1, 2)].r ≈ 1
    @test prob.cexs[CexKey(2, 2)].state.loc == 1
    @test prob.cexs[CexKey(2, 2)].state.x[2] ≈ -1
    @test prob.cexs[CexKey(2, 2)].r ≈ 2
end

nothing