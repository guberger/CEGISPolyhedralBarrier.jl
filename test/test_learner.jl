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
Polyhedron = CPB.Polyhedron
System = CPB.System
InitSet = CPB.InitSet
UnsafeSet = CPB.UnsafeSet

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

## Learner Disc
nvar = 1
nloc = 2

sys = System()

domain = Polyhedron()
CPB.add_halfspace!(domain, [-1], 0)
CPB.add_halfspace!(domain, [1], -2)
A = reshape([-1], 1, 1)
b = [1]
CPB.add_piece!(sys, domain, 1, A, b, 2)

iset = InitSet()
CPB.add_state!(iset, 1, [-1])
CPB.add_state!(iset, 1, [1])

uset = UnsafeSet()
udom = Polyhedron()
CPB.add_halfspace!(udom, [-1], 1.1)
CPB.add_region!(uset, 2, udom)

lear = CPB.Learner(nvar, nloc, sys, iset, uset, 0, 0)
CPB.set_tol!(lear, :rad, 0)
CPB.set_tol!(lear, :bigM, 1e3)

status, = CPB.learn_lyapunov!(lear, 1, solver, solver)

@testset "learn lyapunov disc: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
end

status, = CPB.learn_lyapunov!(lear, 3, solver, solver)

@testset "learn lyapunov disc: found" begin
    @test status == CPB.BARRIER_FOUND
end

domain = Polyhedron()
CPB.add_halfspace!(domain, [1], 0)
A = reshape([-1], 1, 1)
b = [0]
CPB.add_piece!(sys, domain, 1, A, b, 2)

status, = CPB.learn_lyapunov!(lear, 4, solver, solver)

@testset "learn lyapunov disc: found" begin
    @test status == CPB.BARRIER_FOUND
end

CPB.set_tol!(lear, :rad, 0.05)
status, = CPB.learn_lyapunov!(lear, 4, solver, solver)

@testset "learn lyapunov disc: radius too small" begin
    @test status == CPB.RADIUS_TOO_SMALL
end

nothing