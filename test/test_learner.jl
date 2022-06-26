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
Cone = CPB.Cone
System = CPB.System

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

## Learner Cont
nvar = 1
nloc = 1

sys = System()

domain = Cone()
CPB.add_supp!(domain, [1.0])
A = reshape([-1.0], 1, 1)
CPB.add_piece_cont!(sys, domain, 1, A)

domain = Cone()
CPB.add_supp!(domain, [-1.0])
A = reshape([-1.0], 1, 1)
CPB.add_piece_cont!(sys, domain, 1, A)

ϵ = 10.0
δ = 1.0 - 1e-5
τ = 0.1
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 1, solver, solver)[1]

@testset "learn lyapunov cont: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
end

tracerec = CPB.TraceRecorder()
status = CPB.learn_lyapunov!(lear, 2, solver, solver, tracerec=tracerec)[1]

@testset "learn lyapunov cont: found" begin
    @test status == CPB.LYAPUNOV_FOUND
    @test length(tracerec.mpf_list) == 2
    @test length(tracerec.pos_evids_list[2]) == 2
    @test length(tracerec.liedisc_evids_list[2]) == 0
    @test length(tracerec.liecont_evids_list[2]) == 2
end

ϵ = 10.0
δ = 1.0 + 1e-5
τ = 0.1
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 1, solver, solver)[1]

@testset "learn lyapunov cont: infeasible" begin
    @test status == CPB.LYAPUNOV_INFEASIBLE
end

ϵ = 10.0
δ = 1.0 - 1e-5
τ = 0.1
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.1 + 1e-5)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov cont: found" begin
    @test status == CPB.LYAPUNOV_FOUND
end

ϵ = 10.0
δ = 1.0 - 1e-5
τ = 0.1
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.1 + 1e-5)
point = [1.0]
CPB.add_witness!(lear, 1, point)
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov cont: radius too small" begin
    @test status == CPB.RADIUS_TOO_SMALL
end

## Learner Disc
nvar = 1
nloc = 2

sys = System()

domain = Cone()
CPB.add_supp!(domain, [-1.0])
A = reshape([-0.5], 1, 1)
CPB.add_piece_disc!(sys, domain, 1, A, 2)

domain = Cone()
CPB.add_supp!(domain, [1.0])
A = reshape([0.0], 1, 1)
CPB.add_piece_disc!(sys, domain, 2, A, 1)

ϵ = 10.0
δ = min(2/3, 1 - 3/(2*ϵ)) - 1e-5
lear = CPB.Learner(nvar, nloc, sys, 0.0, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 1, solver, solver)[1]

@testset "learn lyapunov disc: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
end

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov disc: found" begin
    @test status == CPB.LYAPUNOV_FOUND
end

ϵ = 10.0
δ = min(2/3, 1 - 3/(2*ϵ)) + 1e-5
lear = CPB.Learner(nvar, nloc, sys, 0.0, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov cont: infeasible" begin
    @test status == CPB.LYAPUNOV_INFEASIBLE
end

ϵ = 10.0
δ = min(2/3, 1 - 3/(2*ϵ)) - 1e-5
lear = CPB.Learner(nvar, nloc, sys, 0.0, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.5 - 1e-5)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov cont: feasible" begin
    @test status == CPB.LYAPUNOV_FOUND
end

ϵ = 10.0
δ = min(2/3, 1 - 3/(2*ϵ)) - 1e-5
lear = CPB.Learner(nvar, nloc, sys, 0.0, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.5 + 1e-5)
point = [1.0]
CPB.add_witness!(lear, 1, point)

status = CPB.learn_lyapunov!(lear, 2, solver, solver)[1]

@testset "learn lyapunov cont: radius too small" begin
    @test status == CPB.RADIUS_TOO_SMALL
end

## Learner Disc and Cont

nvar = 1
nloc = 2

sys = System()

domain = Cone()
CPB.add_supp!(domain, [-1.0])
A = reshape([2.0], 1, 1)
CPB.add_piece_disc!(sys, domain, 1, A, 2)

domain = Cone()
CPB.add_supp!(domain, [-1.0])
A = reshape([-1.0], 1, 1)
CPB.add_piece_cont!(sys, domain, 2, A)

τ = 0.1
ϵ = 10.0
δ = min(1/3, 1 - 1/(2*ϵ)) - 1e-5
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)

status = CPB.learn_lyapunov!(lear, 1, solver, solver)[1]

@testset "learn lyapunov disc & cont: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
end

tracerec = CPB.TraceRecorder()
status = CPB.learn_lyapunov!(lear, 3, solver, solver, tracerec=tracerec)[1]

@testset "learn lyapunov disc & cont: found" begin
    @test status == CPB.LYAPUNOV_FOUND
    @test length(tracerec.mpf_list) == 3
    @test length(tracerec.pos_evids_list[3]) == 2
    @test length(tracerec.liedisc_evids_list[3]) == 1
    @test length(tracerec.liecont_evids_list[3]) == 1
end

τ = 0.1
ϵ = 10.0
δ = min(1/3, 1 - 1/(2*ϵ)) + 1e-5
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 0.0)

status = CPB.learn_lyapunov!(lear, 3, solver, solver)[1]

@testset "learn lyapunov disc & cont: infeasible" begin
    @test status == CPB.LYAPUNOV_INFEASIBLE
end

nothing