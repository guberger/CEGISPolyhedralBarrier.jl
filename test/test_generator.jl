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
NegEvidence = CPB.NegEvidence
PosEvidence = CPB.PosEvidence
LieEvidence = CPB.LieEvidence
PolyFunc = CPB.PolyFunc
_norm(pf::PolyFunc) = maximum(af -> norm(af.a, Inf), pf.afs)

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

βmax = 100

## Empty
gen = CPB.Generator{2}((0,))
M = 1

rf = CPB.compute_mpf_feasibility(gen, 1e-5, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf empty" begin
    @test rf > 0
    @test all(pf -> isempty(pf.afs), mpf.pfs)
    @test re ≈ 10
end

## Pos
gen = CPB.Generator{2}((1,))
M = 8

CPB.add_evidence!(gen, NegEvidence(1, SVector(0.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(1, 1, SVector(0.5, 0.0)))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf pos" begin
    @test rf ≈ 0.5/2 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 0.5/2
end

## Lie
gen = CPB.Generator{2}((1,))
M = 8

CPB.add_evidence!(gen, NegEvidence(1, SVector(0.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(1, 1, SVector(8.0, 0.0)))
point1 = @SVector [2.0, 0.0]
point2 = @SVector [4.0, 0.0]
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 1, point2, 0.0))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf lie" begin
    @test rf ≈ 2 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 4
end

## Lie
gen = CPB.Generator{2}((1,))
M = 32

CPB.add_evidence!(gen, NegEvidence(1, SVector(0.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(1, 1, SVector(8.0, 0.0)))
point1 = @SVector [6.0, 0.0]
point2 = @SVector [4.0, 0.0]
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 1, point2, 0.0))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf lie" begin
    @test rf ≈ 3 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 4
end

## Pos and Lie: 2 wits #1
gen = CPB.Generator{2}((1, 1))
M = 32

CPB.add_evidence!(gen, NegEvidence(1, SVector(-2.0, 0.0)))
CPB.add_evidence!(gen, NegEvidence(1, SVector(2.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(2, 1, SVector(2.0, 0.0)))
point1 = @SVector [2.0, 0.0]
point2 = @SVector [4.0, 0.0]
nA = 0.5
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 2, point2, nA))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf pos and lie: 2 wits #1" begin
    @test rf ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 2/(1 + nA)
end

## Pos and Lie: 2 wits #2
gen = CPB.Generator{2}((1, 1))
M = 32

CPB.add_evidence!(gen, NegEvidence(1, SVector(-2.0, 0.0)))
CPB.add_evidence!(gen, NegEvidence(1, SVector(2.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(2, 1, SVector(2.0, 0.0)))
point1 = @SVector [4.0, 0.0]
point2 = @SVector [2.0, 0.0]
nA = 0.5
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 2, point2, nA))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf pos and lie: 2 wits #2" begin
    @test rf ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 1
end

## Pos and Lie: 2 wits #3
gen = CPB.Generator{2}((2, 1))
M = 32

CPB.add_evidence!(gen, NegEvidence(1, SVector(-2.0, 0.0)))
CPB.add_evidence!(gen, NegEvidence(1, SVector(2.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(1, 1, SVector(-4.0, 0.0)))
CPB.add_evidence!(gen, PosEvidence(1, 2, SVector(4.0, 0.0)))
point1 = @SVector [1.0, 0.0]
point2 = @SVector [1.0, 0.0]
point3 = @SVector [2.0, 0.0]
nA = 1.0
CPB.add_evidence!(gen, LieEvidence(1, 2, point1, 2, point2, nA))
CPB.add_evidence!(gen, LieEvidence(2, 1, point2, 1, point3, nA))

ϵ = 1e-2
rf = CPB.compute_mpf_feasibility(gen, ϵ, M, βmax, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, βmax, solver)

@testset "compute mpf pos and lie: 2 wits #3" begin
    @test rf ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 1
end

nothing