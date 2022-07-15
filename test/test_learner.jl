using LinearAlgebra
using StaticArrays
using JuMP
using HiGHS
using Test
@static if isdefined(Main, :TestLocal)
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier
System = CPB.System
PointSet = CPB.PointSet
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
_norm(pf::PolyFunc) = maximum(af -> norm(af.a, Inf), pf.afs)

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

sys = System{1}()

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(1.0), -2.0)
A = @SMatrix [0.5]
b = @SVector [0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 2)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(1.0), -2.0)
A = @SMatrix [1.0]
b = @SVector [0.5]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

iset = PointSet{1,2}()
CPB.add_point!(iset, 1, SVector(0.5))

mpf_safe = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_safe, 1, SVector(1.0), -1.8)

mpf_inv = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0), 0.0)
CPB.add_af!(mpf_inv, 2, SVector(-1.0), 0.0)

lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, 1e-3, 1e-8)
CPB.set_param!(lear, :xmax, 1e2)

@testset "set tol and param" begin
    @test_throws AssertionError CPB.set_param!(lear, :dumb, 0)
    @test lear.params[:xmax] ≈ 100
end

lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, 1e-1, 1e-8)
status, = CPB.learn_lyapunov!(lear, 1, solver, solver)

@testset "learn lyapunov disc: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
end

status, = CPB.learn_lyapunov!(
    lear, 500, solver, solver, do_print=true
)

@testset "learn lyapunov disc: found" begin
    @test status == CPB.BARRIER_FOUND
end

lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, 0.55, 1e-8)
status, = CPB.learn_lyapunov!(lear, 20, solver, solver)

@testset "learn lyapunov disc: radius too small" begin
    @test status == CPB.RADIUS_TOO_SMALL
end

nothing