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

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

xmax = 1e3
η = -0.1

# Set #1
N = 1
pf_dom = PolyFunc([AffForm([-1.0], 0.0)])
A = [0.5;;]
b = [1.0]
piece = Piece(pf_dom, 1, A, b, 1)
sys = System([piece])

mpf_inv = MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])])

mpf_BF = MultiPolyFunc([PolyFunc([AffForm([-0.5], 0.25)])])

mpf_safe = MultiPolyFunc([PolyFunc(AffForm{Vector{Int},Int}[])])

x, r, loc = CPB.verify_safe(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify safe: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

mpf_safe = MultiPolyFunc([PolyFunc([AffForm([1.0], -0.25)])])

x, r, loc = CPB.verify_safe(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify safe: infeasible" begin
    @test r ≈ -0.125
    @test x ≈ [0.375]
    @test loc == 1
end

mpf_safe = MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])])

x, r, loc = CPB.verify_safe(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify safe: unsafe" begin
    @test r ≈ 0.25
    @test x ≈ [0.75]
    @test loc == 1
end

x, r, loc = CPB.verify_BF(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify BF: satisfied" begin
    @test r ≈ -2*η/3 - 0.5
    @test x ≈ [-2*η/3]
    @test loc == 1
end

# Set #2
N = 2
pf_dom = PolyFunc([AffForm([-1.0, 0.0], -1.0)])
A = [0.0 1.0; -1.0 0.0]
b = [0.0, 0.0]
sys = System([Piece(pf_dom, 1, A, b, 2)])

mpf_inv = MultiPolyFunc([
    PolyFunc([AffForm([0.0, 1.0], -1.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

mpf_BF = MultiPolyFunc([
    PolyFunc([AffForm([1.0, 0.0], -1.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

mpf_safe = MultiPolyFunc([
    PolyFunc([AffForm([0.0, -1.0], -1.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

x, r, loc = CPB.verify_safe(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify safe: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

x, r, loc = CPB.verify_BF(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify BF: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

mpf_BF = MultiPolyFunc([
    PolyFunc(AffForm{Vector{Float64},Float64}[]),
    PolyFunc([AffForm([1.0, 0.0], -1.0)])
])

x, r, loc = CPB.verify_BF(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify BF: unsatisfied" begin
    @test r ≈ -η
    @test x[2] ≈ 1
    @test loc == 1
end

pf_dom = PolyFunc([AffForm([-1.0, 0.0], -1.0)])
A = [0.0 1.0; -1.0 0.0]
b = [0.0, 0.0]
piece1 = Piece(pf_dom, 1, A, b, 2)
pf_dom = PolyFunc([AffForm([-1.0, 0.0], -1.0)])
A = [0.0 1.0; -1.0 0.0]
b = [0.0, -1.0]
piece2 = Piece(pf_dom, 1, A, b, 2)
sys = System([piece1, piece2])

mpf_inv = MultiPolyFunc([
    PolyFunc([AffForm([1.0, 0.0], -2.0)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

mpf_BF = MultiPolyFunc([
    PolyFunc(AffForm{Vector{Float64},Float64}[]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

mpf_safe = MultiPolyFunc([
    PolyFunc([AffForm([0.0, -1.0], -1.0)]),
    PolyFunc([AffForm([0.0, -1.0], -1.0)])
])

x, r, loc = CPB.verify_safe(sys, mpf_safe, mpf_inv, mpf_BF, xmax, η, N, solver)

@testset "verify safe: unsafe" begin
    @test r ≈ xmax + 1
    @test x[1] ≈ 1 + η
    @test loc == 1
end

nothing