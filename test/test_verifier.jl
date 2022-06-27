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
AffForm = CPB.AffForm
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
PosPredicate = CPB.PosPredicate
LiePredicate = CPB.LiePredicate

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

## Parameters
nvar = 2

## Pos false #1 #2 #3
verif = CPB.Verifier()
domain = Polyhedron()
CPB.add_halfspace!(domain, [1.0, 1.0], -1)
CPB.add_predicate!(verif, PosPredicate(nvar, domain, 1))

mpf = MultiPolyFunc(2)
afs = [AffForm([-0.5, 0.5], -0.5), AffForm([1.0, 0.0], 0)]
for af in afs
    CPB.add_af!(mpf, 1, af)
end

xmax = 1e3
r, x, loc = CPB.verify_pos(verif, mpf, xmax, 1e4, solver)

@testset "verify pos false #1" begin
    @test r ≈ (xmax + 1)/3
    @test norm(x, Inf) ≈ xmax
    @test x ∈ domain
    @test loc == 1
end

domain = Polyhedron()
CPB.add_af!(mpf, 2, AffForm([4.0, 0.0], 1))
CPB.add_predicate!(verif, PosPredicate(nvar, domain, 2))

r, x, loc = CPB.verify_pos(verif, mpf, xmax, 1e4, solver)

@testset "verify pos false #2" begin
    @test r ≈ 4*xmax - 1
    @test loc == 2
end

CPB.add_af!(mpf, 2, AffForm([4.0, 0.0], 4000))

r, x, loc = CPB.verify_pos(verif, mpf, xmax, 1e4, solver)

@testset "verify pos false #3" begin
    @test r ≈ (xmax + 1)/3
    @test loc == 1
end

# Pos true
verif = CPB.Verifier()
domain = Polyhedron()
CPB.add_halfspace!(domain, [-1.0, -1.0], 1)
CPB.add_predicate!(verif, PosPredicate(nvar, domain, 1))

mpf = MultiPolyFunc(2)
afs = [AffForm([-0.5, 0.5], 0.5), AffForm([1.0, 0.0], 0)]
for af in afs
    CPB.add_af!(mpf, 1, af)
end

xmax = 1e3
r, x, loc = CPB.verify_pos(verif, mpf, xmax, 1e4, solver)

@testset "verify pos true" begin
    @test r ≈ -1/2
    @test x ≈ [1/2, 1/2]
    @test x ∈ domain
    @test loc == 1
end

## Lie disc false #1
verif = CPB.Verifier()
domain = Polyhedron()
A = [0.5 0.0; 1.0 1.0]
b = [1, 0]
CPB.add_predicate!(verif, LiePredicate(nvar, domain, 1, A, b, 1))

mpf = MultiPolyFunc(2)
afs = [AffForm([-1.0, 0.0], -1), AffForm([1.0, 0.0], -1)]
for af in afs
    CPB.add_af!(mpf, 1, af)
end

r, x, loc = CPB.verify_lie(verif, mpf, 1e3, 1e4, solver)

@testset "verify lie disc false #1" begin
    @test r ≈ 1/2
    @test x[1] ≈ 1
    @test x ∈ domain
    @test loc == 1
end

## Lie disc false #2
verif = CPB.Verifier()
domain = Polyhedron()
CPB.add_halfspace!(domain, [-1.0, 0.0], 0)
A = [0.0 0.5; 0.5 0.1]
b = [0, 1]
CPB.add_predicate!(verif, LiePredicate(nvar, domain, 1, A, b, 1))

mpf = MultiPolyFunc(2)
afs = [
    AffForm([-1.0, 0.0], -1), AffForm([1.0, 0.0], -1),
    AffForm([0.0, -1.0], -1), AffForm([0.0, 1.0], -1)
]
for af in afs
    CPB.add_af!(mpf, 1, af)
end

r, x, loc = CPB.verify_lie(verif, mpf, 1e3, 1e4, solver)

@testset "verify lie disc false #2" begin
    @test r ≈ 0.6
    @test x ≈ [1, 1]
    @test x ∈ domain
    @test loc == 1
end

## Lie disc true #1
verif = CPB.Verifier()
domain = Polyhedron()
CPB.add_halfspace!(domain, [-1.0, 0.0], 0)
A = [0.0 0.5; 0.5 0.1]
b = [0, -0.5]
CPB.add_predicate!(verif, LiePredicate(nvar, domain, 1, A, b, 1))

mpf = MultiPolyFunc(2)
afs = [
    AffForm([-1.0, 0.0], -1), AffForm([1.0, 0.0], -1),
    AffForm([0.0, -1.0], -1), AffForm([0.0, 1.0], -1)
]
for af in afs
    CPB.add_af!(mpf, 1, af)
end

r, x, loc = CPB.verify_lie(verif, mpf, 1e3, 1e4, solver)

@testset "verify lie disc false #2" begin
    @test r ≈ -0.4
    @test x ≈ [0, -1]
    @test x ∈ domain
    @test loc == 1
end

## Lie disc multiple #1
verif = CPB.Verifier()
domain = Polyhedron()
CPB.add_halfspace!(domain, [-1.0, -1.0], -1)
CPB.add_halfspace!(domain, [-1.0, 1.0], -1)
A = [0.5 -0.25; 0.1 0.5]
b = [0, 0.5]
CPB.add_predicate!(verif, LiePredicate(nvar, domain, 2, A, b, 1))

mpf = MultiPolyFunc(2)
afs = [
    AffForm([1.0, 0.0], -1), AffForm([0.0, -1.0], -1), AffForm([0.0, 1.0], -1)
]
for af in afs
    CPB.add_af!(mpf, 1, af)
    CPB.add_af!(mpf, 2, af)
end

r, x, loc = CPB.verify_lie(verif, mpf, 1e3, 1e4, solver)

@testset "verify lie disc multiple #1" begin
    @test r ≈ 0.1
    @test x ≈ [1, 1]
    @test x ∈ domain
    @test loc == 2
end

mpf = MultiPolyFunc(2)
afs = [
    AffForm([1.0, 0.0], -1), AffForm([0.0, -1.0], -1), AffForm([0.0, 1.0], -1)
]
for af in afs
    CPB.add_af!(mpf, 2, af)
end
CPB.add_af!(mpf, 1, afs[1])

r, x, loc = CPB.verify_lie(verif, mpf, 1e3, 1e4, solver)

@testset "verify lie disc multiple #2" begin
    @test r ≈ -0.25
    @test x ≈ [1, -1]
    @test x ∈ domain
    @test loc == 2
end

nothing