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
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
Piece = CPB.Piece
System = CPB.System
Witness = CPB.Witness

function HiGHS._check_ret(ret::Cint) 
    if ret != Cint(0) && ret != Cint(1)
        error(
            "Encountered an error in HiGHS (Status $(ret)). Check the log " * 
            "for details.", 
        )
    end 
    return 
end 

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

#-------------------------------------------------------------------------------
N = 1
M = 1
pf_dom = PolyFunc([AffForm([-1.0], 3.0), AffForm([1.0], -3.0)])
A = [1.0;;]
b = [1.0]
sys = System([Piece(pf_dom, 1, A, b, 1)])

mlist_init = [[[1.0]]]

mpf_safe = MultiPolyFunc([PolyFunc([AffForm([1.0], -4.0)])])
mpf_inv = MultiPolyFunc([PolyFunc([AffForm([-1.0], 0.0)])])

wit_trace = Witness[]
function rec_wit_trace(::Any, ::Any, wit)
    push!(wit_trace, Witness(
        copy.(wit.mlist_inside), copy.(wit.mlist_image),
        copy.(wit.mlist_unknown), copy.(wit.mlist_outside)
    ))
end

ϵ = 100.0
δ = 1e-3
xmax = 1e2
tol_dom = 1e-6
iter_max = 1

status, = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
    @test length(wit_trace) == 1
end

iter_max = 3
empty!(wit_trace)

status, mpf, wit_final = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: feasible" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wit_trace) == 2
    #1
    @test wit_trace[1].mlist_inside[1] ≈ [[1]]
    @test isempty(wit_trace[1].mlist_image[1])
    @test isempty(wit_trace[1].mlist_outside[1])
    @test isempty(wit_trace[1].mlist_unknown[1])
    #2
    @test wit_trace[2].mlist_inside[1] ≈ [[1]]
    @test isempty(wit_trace[2].mlist_image[1])
    @test wit_trace[2].mlist_outside[1] ≈ [[3]]
    @test isempty(wit_trace[2].mlist_unknown[1])
end

pf_dom = PolyFunc([AffForm([-1.0], 3.0), AffForm([1.0], -3.0)])
A = [1.0;;]
b = [1.0]
piece1 = Piece(pf_dom, 1, A, b, 1)
pf_dom = PolyFunc([AffForm([-1.0], 1.0), AffForm([1.0], -1.0)])
A = [0.0;;]
b = [3.5]
piece2 = Piece(pf_dom, 1, A, b, 1)
sys = System([piece1, piece2])

ϵ = 1e-4
δ = 1e-3
xmax = 1e2
tol_dom = 1e-6
iter_max = 3
empty!(wit_trace)

status, mpf, wit_final = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: rad too small" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(wit_trace) == 2
    #1
    @test wit_trace[1].mlist_inside[1] ≈ [[1]]
    @test wit_trace[1].mlist_image[1] ≈ [[3.5]]
    @test isempty(wit_trace[1].mlist_outside[1])
    @test isempty(wit_trace[1].mlist_unknown[1])
    #2
    @test wit_trace[2].mlist_inside[1] ≈ [[1]]
    @test wit_trace[2].mlist_image[1] ≈ [[3.5]]
    @test wit_trace[2].mlist_outside[1] ≈ [[3]]
    @test isempty(wit_trace[2].mlist_unknown[1])
end

#-------------------------------------------------------------------------------
N = 1
M = 2
pf_dom = PolyFunc([AffForm([-1.0], 3.0), AffForm([1.0], -3.0)])
A = [0.5;;]
b = [2.5]
piece1 = Piece(pf_dom, 1, A, b, 1)
pf_dom = PolyFunc([AffForm([-1.0], 1.0), AffForm([1.0], -1.0)])
A = [0.0;;]
b = [2.5]
piece2 = Piece(pf_dom, 2, A, b, 1)
pf_dom = PolyFunc([AffForm([1.0], -0.0)])
A = [0.0;;]
b = [0.9]
piece3 = Piece(pf_dom, 2, A, b, 2)
sys = System([piece1, piece2, piece3])

mlist_init = [[[1.0]], [[0.0]]]

mpf_safe = MultiPolyFunc([
    PolyFunc([AffForm([1.0], -4.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])
mpf_inv = MultiPolyFunc([
    PolyFunc([AffForm([-1.0], 0.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

ϵ = 100.0
δ = 1e-3
xmax = 1e2
tol_dom = 1e-6
iter_max = 4
empty!(wit_trace)

status, mpf, wit_final = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: unknown -> outside" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(wit_trace) == 4
    #1
    @test wit_trace[4].mlist_inside[1] ≈ [[1]]
    @test wit_trace[4].mlist_inside[2] ≈ [[0]]
    @test isempty(wit_trace[4].mlist_image[1])
    @test wit_trace[4].mlist_image[2] ≈ [[0.9]]
    @test wit_trace[4].mlist_outside[1] ≈ [[3]]
    @test wit_trace[4].mlist_outside[2] ≈ [[1]]
    @test isempty(wit_trace[4].mlist_unknown[1])
    @test isempty(wit_trace[4].mlist_unknown[2])
end

ϵ = 1e-2
δ = 1e-3
xmax = 1e2
tol_dom = 1e-6
iter_max = 4
empty!(wit_trace)

status, mpf, wit_final = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: unknown! gets excluded" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wit_trace) == 4
    #1
    @test wit_trace[4].mlist_inside[1] ≈ [[1]]
    @test wit_trace[4].mlist_inside[2] ≈ [[0]]
    @test isempty(wit_trace[4].mlist_image[1])
    @test wit_trace[4].mlist_image[2] ≈ [[0.9]]
    @test wit_trace[4].mlist_outside[1] ≈ [[3]]
    @test isempty(wit_trace[4].mlist_outside[2])
    @test isempty(wit_trace[4].mlist_unknown[1])
    @test wit_trace[4].mlist_unknown[2] ≈ [[1]]
end

ϵ = 1e-1
δ = 1e-3
xmax = 1e2
tol_dom = 1e-6
iter_max = 6
empty!(wit_trace)

status, mpf, wit_final = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver,
    tol_dom=tol_dom, xmax=xmax, callback_fcn=rec_wit_trace
)

@testset "learn lyapunov: unknown! gets included" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wit_trace) == 6
    #1
    @test wit_trace[6].mlist_inside[1] ≈ [[1]]
    @test wit_trace[6].mlist_inside[2] ≈ [[0], [1]]
    @test wit_trace[6].mlist_image[1] ≈ [[2.5]]
    @test wit_trace[6].mlist_image[2] ≈ [[0.9]]
    @test wit_trace[6].mlist_outside[1] ≈ [[3]]
    @test isempty(wit_trace[6].mlist_outside[2])
    @test isempty(wit_trace[6].mlist_unknown[1])
    @test isempty(wit_trace[6].mlist_unknown[2])
end

nothing