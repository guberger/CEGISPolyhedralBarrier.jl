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

# Set 1
sys = System{1}()
ϵ = 0.5 - 1e-3

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(1.0), -2.0)
A = @SMatrix [0.5]
b = @SVector [0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 2.0)
CPB.add_af!(pf_dom, SVector(1.0), -3.0)
A = @SMatrix [0.0]
b = @SVector [4 - ϵ/2]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 3.0)
A = @SMatrix [1.0]
b = @SVector [ϵ/32]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

iset = PointSet{1,2}()
CPB.add_point!(iset, 1, SVector(1.0))

mpf_safe = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_safe, 1, SVector(1.0), -4 - ϵ/16)

mpf_inv = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0), 0.0)

lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, ϵ, ϵ/32)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, ϵ/1024)

@testset "set tol and param" begin
    @test_throws AssertionError CPB.set_param!(lear, :dumb, 0)
    @test lear.params[:xmax] ≈ 100
end

wits = CPB.Witness{1,2}[]
function rec_wits(::Any, ::Any, wit)
    inside_ = PointSet(copy.(wit.inside.points_list))
    image_ = PointSet(copy.(wit.image.points_list))
    outside_ = PointSet(copy.(wit.outside.points_list))
    unknown_ = PointSet(copy.(wit.unknown.points_list))
    push!(wits, CPB.Witness(inside_, image_, outside_, unknown_))
end
status, = CPB.learn_lyapunov!(lear, 1, solver, solver, callback_fcn=rec_wits)

@testset "learn lyapunov #1: max iter" begin
    @test status == CPB.MAX_ITER_REACHED
    @test length(wits) == 1
end

wits = CPB.Witness{1,2}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 3, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov #2: found" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wits) == 3
    #1
    @test wits[1].inside.points_list[1] ≈ [[1]]
    @test isempty(wits[1].inside.points_list[2])
    @test wits[1].image.points_list[1] ≈ [[0.5]]
    @test isempty(wits[1].image.points_list[2])
    @test isempty(wits[1].outside.points_list[1])
    @test isempty(wits[1].outside.points_list[2])
    @test isempty(wits[1].unknown.points_list[1])
    @test isempty(wits[1].unknown.points_list[2])
    #2
    @test wits[2].inside.points_list[1] ≈ [[1]]
    @test isempty(wits[2].inside.points_list[2])
    @test wits[2].image.points_list[1] ≈ [[0.5]]
    @test isempty(wits[2].image.points_list[2])
    @test wits[2].outside.points_list[1] ≈ [[4]]
    @test isempty(wits[2].outside.points_list[2])
    @test isempty(wits[2].unknown.points_list[1])
    @test isempty(wits[2].unknown.points_list[2])
    #3
    @test wits[3].inside.points_list[1] ≈ [[1]]
    @test isempty(wits[3].inside.points_list[2])
    @test wits[3].image.points_list[1] ≈ [[0.5]]
    @test isempty(wits[3].image.points_list[2])
    @test wits[3].outside.points_list[1] ≈ [[4], [2.125]]
    @test isempty(wits[3].outside.points_list[2])
    @test isempty(wits[3].unknown.points_list[1])
    @test isempty(wits[3].unknown.points_list[2])
    #mpf
    @test any(af -> af.a ≈ [1] && af.β ≈ -1.3125, mpf.pfs[1].afs)
    @test isempty(mpf.pfs[2].afs)
end

# Set 2
sys = System{1}()

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(1.0), 0.0)
A = @SMatrix [1.0]
b = @SVector [0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 2)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 1.0)
A = @SMatrix [1.0]
b = @SVector [1.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 0.25)
CPB.add_af!(pf_dom, SVector(1.0), -0.25)
A = @SMatrix [0.0]
b = @SVector [0.75]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

iset = PointSet{1,2}()
CPB.add_point!(iset, 1, SVector(0.0))

mpf_safe = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_safe, 1, SVector(1.0), -2.0)

mpf_inv = MultiPolyFunc{1,2}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0), 0.0)
CPB.add_af!(mpf_inv, 2, SVector(-1.0), 0.0)

lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, 0.125 + 1e-3, 1e-8)
CPB.set_param!(lear, :xmax, 1e2)
CPB.set_param!(lear, :tol_dom, 1e-7)
wits = CPB.Witness{1,2}[]
function rec_wits(::Any, ::Any, wit)
    inside_ = PointSet(copy.(wit.inside.points_list))
    image_ = PointSet(copy.(wit.image.points_list))
    outside_ = PointSet(copy.(wit.outside.points_list))
    unknown_ = PointSet(copy.(wit.unknown.points_list))
    push!(wits, CPB.Witness(inside_, image_, outside_, unknown_))
end
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 7, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov #2: found" begin
    @test status == CPB.BARRIER_FOUND
    @test length(wits) == 7
    #1
    @test wits[1].inside.points_list[1] ≈ [[0]]
    @test isempty(wits[1].inside.points_list[2])
    @test isempty(wits[1].image.points_list[1])
    @test wits[1].image.points_list[2] ≈ [[0.0]]
    @test isempty(wits[1].outside.points_list[1])
    @test isempty(wits[1].outside.points_list[2])
    @test isempty(wits[1].unknown.points_list[1])
    @test isempty(wits[1].unknown.points_list[2])
    #2
    @test wits[2].inside.points_list[1] ≈ [[0]]
    @test isempty(wits[2].inside.points_list[2])
    @test isempty(wits[2].image.points_list[1])
    @test wits[2].image.points_list[2] ≈ [[0.0]]
    @test isempty(wits[2].outside.points_list[1])
    @test isempty(wits[2].outside.points_list[2])
    @test isempty(wits[2].unknown.points_list[1])
    @test isempty(wits[2].unknown.points_list[2])
    #3
    @test wits[3].inside.points_list[1] ≈ [[0]]
    @test isempty(wits[3].inside.points_list[2])
    @test isempty(wits[3].image.points_list[1])
    @test wits[3].image.points_list[2] ≈ [[0.0]]
    @test wits[3].outside.points_list[1] ≈ [[1.5]]
    @test isempty(wits[3].outside.points_list[2])
    @test isempty(wits[3].unknown.points_list[1])
    @test isempty(wits[3].unknown.points_list[2])
    #4
    @test wits[4].inside.points_list[1] ≈ [[0]]
    @test isempty(wits[4].inside.points_list[2])
    @test isempty(wits[4].image.points_list[1])
    @test wits[4].image.points_list[2] ≈ [[0.0]]
    @test wits[4].outside.points_list[1] ≈ [[1.5]]
    @test isempty(wits[4].outside.points_list[2])
    @test isempty(wits[4].unknown.points_list[1])
    @test isapprox(wits[4].unknown.points_list[2], [[0.25]], atol=1e-6)
    #mpf
    @test any(af -> af.a ≈ [1] && af.β ≈ -0.90625, mpf.pfs[1].afs)
    @test isempty(mpf.pfs[2].afs)
end

pf_dom = PolyFunc{1}()
CPB.add_af!(pf_dom, SVector(-1.0), 0.25)
CPB.add_af!(pf_dom, SVector(1.0), -0.25)
A = @SMatrix [0.0]
b = @SVector [2.0]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

wits = CPB.Witness{1,2}[]
status, mpf, wit = CPB.learn_lyapunov!(
    lear, 7, solver, solver, callback_fcn=rec_wits
)

@testset "learn lyapunov #2: infeasible" begin
    @test status == CPB.BARRIER_INFEASIBLE
    @test length(wits) == 4
end

nothing