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
MultiSet = CPB.MultiSet
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

#-------------------------------------------------------------------------------
sys = System{1}()

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 3.0)
CPB.add_af!(pf_dom, SVector(1.0), -3.0)
A = @SMatrix [1.0]
b = @SVector [1.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

mset_init = MultiSet{1,1}()
CPB.add_point!(mset_init, 1, SVector(1.0))

mpf_safe = MultiPolyFunc{1,1}()
CPB.add_af!(mpf_safe, 1, SVector(1.0), -4.0)

mpf_inv = MultiPolyFunc{1,1}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0), 0.0)

ϵ = 100.0
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, ϵ, 1e-3)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-6)

@testset "set tol and param" begin
    @test_throws AssertionError CPB.set_param!(lear, :dumb, 0)
    @test lear.params[:xmax] ≈ 100
    @test lear.params[:tol_dom] ≈ 1e-6
end

wits = CPB.Witness{1,1}[]
function rec_wits(::Any, ::Any, wit)
    inside_ = MultiSet(copy.(wit.inside.sets))
    image_ = MultiSet(copy.(wit.image.sets))
    outside_ = MultiSet(copy.(wit.outside.sets))
    unknown_ = MultiSet(copy.(wit.unknown.sets))
    push!(wits, CPB.Witness(inside_, image_, outside_, unknown_))
end
status, = CPB.learn_lyapunov!(lear, 1, solver, solver, callback_fcn=rec_wits)

@testset "learn lyapunov: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
    @test length(wits) == 1
end

wits = CPB.Witness{1,1}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 3, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov: feasible" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wits) == 2
    #1
    @test wits[1].inside.sets[1] ≈ [[1]]
    @test isempty(wits[1].image.sets[1])
    @test isempty(wits[1].outside.sets[1])
    @test isempty(wits[1].unknown.sets[1])
    #2
    @test wits[2].inside.sets[1] ≈ [[1]]
    @test isempty(wits[2].image.sets[1])
    @test wits[2].outside.sets[1] ≈ [[3]]
    @test isempty(wits[2].unknown.sets[1])
end

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 1.0)
CPB.add_af!(pf_dom, SVector(1.0), -1.0)
A = @SMatrix [0.0]
b = @SVector [3.5]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

ϵ = 1e-4
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, ϵ, 1e-3)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-6)

wits = CPB.Witness{1,1}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 3, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov: rad too small" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(wits) == 2
    #1
    @test wits[1].inside.sets[1] ≈ [[1]]
    @test wits[1].image.sets[1] ≈ [[3.5]]
    @test isempty(wits[1].outside.sets[1])
    @test isempty(wits[1].unknown.sets[1])
    #2
    @test wits[2].inside.sets[1] ≈ [[1]]
    @test wits[2].image.sets[1] ≈ [[3.5]]
    @test wits[2].outside.sets[1] ≈ [[3]]
    @test isempty(wits[2].unknown.sets[1])
end

#-------------------------------------------------------------------------------
sys = System{1}()

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 3.0)
CPB.add_af!(pf_dom, SVector(1.0), -3.0)
A = @SMatrix [0.5]
b = @SVector [2.5]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 1.0)
CPB.add_af!(pf_dom, SVector(1.0), -1.0)
A = @SMatrix [0.0]
b = @SVector [2.5]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(1.0), -0.0)
A = @SMatrix [0.0]
b = @SVector [0.9]
CPB.add_piece!(sys, pf_dom, 2, A, b, 2)

mset_init = MultiSet{1,2}()
CPB.add_point!(mset_init, 1, SVector(1.0))
CPB.add_point!(mset_init, 2, SVector(0.0))

mpf_safe = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_safe, 1, SVector(1.0), -4.0)

mpf_inv = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0), 0.0)

ϵ = 100.0
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, ϵ, 1e-3)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-6)

wits = CPB.Witness{1,2}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 4, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov: unknown -> outside" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(wits) == 4
    #1
    @test wits[4].inside.sets[1] ≈ [[1]]
    @test wits[4].inside.sets[2] ≈ [[0]]
    @test isempty(wits[4].image.sets[1])
    @test wits[4].image.sets[2] ≈ [[0.9]]
    @test wits[4].outside.sets[1] ≈ [[3]]
    @test wits[4].outside.sets[2] ≈ [[1]]
    @test isempty(wits[4].unknown.sets[1])
    @test isempty(wits[4].unknown.sets[2])
end

ϵ = 0.01
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, ϵ, 1e-3)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-6)

wits = CPB.Witness{1,2}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 4, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov: unknown! gets excluded" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wits) == 4
    #1
    @test wits[4].inside.sets[1] ≈ [[1]]
    @test wits[4].inside.sets[2] ≈ [[0]]
    @test isempty(wits[4].image.sets[1])
    @test wits[4].image.sets[2] ≈ [[0.9]]
    @test wits[4].outside.sets[1] ≈ [[3]]
    @test isempty(wits[4].outside.sets[2])
    @test isempty(wits[4].unknown.sets[1])
    @test wits[4].unknown.sets[2] ≈ [[1]]
end

ϵ = 0.1
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, ϵ, 1e-3)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-6)

wits = CPB.Witness{1,2}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 6, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov: unknown! gets included" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wits) == 6
    #1
    @test wits[6].inside.sets[1] ≈ [[1]]
    @test wits[6].inside.sets[2] ≈ [[0], [1]]
    @test wits[6].image.sets[1] ≈ [[2.5]]
    @test wits[6].image.sets[2] ≈ [[0.9]]
    @test wits[6].outside.sets[1] ≈ [[3]]
    @test isempty(wits[6].outside.sets[2])
    @test isempty(wits[6].unknown.sets[1])
    @test isempty(wits[6].unknown.sets[2])
end

nothing