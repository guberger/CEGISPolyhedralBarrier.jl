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
VerifierSafeProblem = CPB.VerifierSafeProblem
VerifierBFProblem = CPB.VerifierBFProblem
empty_pf = CPB.empty_pf
empty_mpf = CPB.empty_mpf

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

xmax = 1e3

# Set #1
N = 1
pf_dom = PolyFunc([AffForm([-1.0], 0.0)])
A = [0.5;;]
b = [1.0]
piece = Piece(pf_dom, 1, A, b, 1)
sys = System([piece])

prob = VerifierSafeProblem(N, sys,
                           empty_mpf(1),
                           MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                           MultiPolyFunc([PolyFunc([AffForm([-0.5], 0.25)])]))

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify safe: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

prob = VerifierSafeProblem(N, sys,
                           MultiPolyFunc([PolyFunc([AffForm([1.0], -0.25)])]),
                           MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                           MultiPolyFunc([PolyFunc([AffForm([-0.5], 0.25)])]))

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify safe: infeasible" begin
    @test r ≈ -0.125
    @test x ≈ [0.375]
    @test loc == 1
end

prob = VerifierSafeProblem(N, sys,
                           MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                           MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                           MultiPolyFunc([PolyFunc([AffForm([-0.5], 0.25)])]))

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify safe: unsafe" begin
    @test r ≈ 0.25
    @test x ≈ [0.75]
    @test loc == 1
end

prob = VerifierBFProblem(N, sys,
                         MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                         MultiPolyFunc([PolyFunc([AffForm([1.0], -1.0)])]),
                         MultiPolyFunc([PolyFunc([AffForm([-0.5], 0.25)])]))

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify BF: satisfied" begin
    @test r ≈ -0.5
    @test norm(x) < 1e-6
    @test loc == 1
end

# Set #2
N = 2
pf_dom = PolyFunc([AffForm([-1.0, 0.0], -1.0)])
A = [0.0 1.0; -1.0 0.0]
b = [0.0, 0.0]
sys = System([Piece(pf_dom, 1, A, b, 2)])

prob = VerifierSafeProblem(
    N, sys,
    MultiPolyFunc([PolyFunc([AffForm([0.0, -1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([PolyFunc([AffForm([0.0, 1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([PolyFunc([AffForm([1.0, 0.0], -1.0)]), empty_pf()])
)

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify safe: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

prob = VerifierBFProblem(
    N, sys,
    MultiPolyFunc([PolyFunc([AffForm([0.0, -1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([PolyFunc([AffForm([0.0, 1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([PolyFunc([AffForm([1.0, 0.0], -1.0)]), empty_pf()])
)

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify BF: empty" begin
    @test r == -Inf
    @test all(isnan, x)
    @test loc == 0
end

prob = VerifierBFProblem(
    N, sys,
    MultiPolyFunc([PolyFunc([AffForm([0.0, -1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([PolyFunc([AffForm([0.0, 1.0], -1.0)]), empty_pf()]),
    MultiPolyFunc([empty_pf(), PolyFunc([AffForm([1.0, 0.0], -1.0)])])
)

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify BF: unsatisfied" begin
    @test abs(r) < 1e-6
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

prob = VerifierSafeProblem(
    N, sys,
    MultiPolyFunc([PolyFunc([AffForm([0.0, -1.0], -1.0)]),
                   PolyFunc([AffForm([0.0, -1.0], -1.0)])]),
    MultiPolyFunc([PolyFunc([AffForm([1.0, 0.0], -2.0)]), empty_pf()]),
    empty_mpf(2)
)

x, r, loc = CPB.find_counterexample(prob, xmax, solver)

@testset "verify safe: unsafe" begin
    @test r ≈ xmax + 1
    @test x[1] ≈ 1
    @test loc == 1
end

nothing