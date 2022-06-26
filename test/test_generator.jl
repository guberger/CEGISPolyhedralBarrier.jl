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
NegEvidence = CPB.NegEvidence
PosEvidence = CPB.PosEvidence
LieEvidence = CPB.LieEvidence
PolyFunc = CPB.PolyFunc
_norm(pf::PolyFunc) = maximum(lf -> norm(lf.lin, 1), pf.afs)

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

## Parameters
nvar = 2

## Empty
nloc = 1
gen = CPB.Generator(nvar, nloc)
M = 1.0

rf = CPB.compute_mpf_feasibility(gen, 1e5, 1.0, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf empty" begin
    @test rf > 0
    @test all(pf -> isempty(pf.afs), mpf.pfs)
    @test re ≈ 1e4
end

## Pos
nloc = 1
gen = CPB.Generator(nvar, nloc)
M = 8.0

CPB.add_af!(gen, 1)
CPB.add_evidence!(gen, NegEvidence(1, [0, 0]))
point = [0.5, 0]
CPB.add_evidence!(gen, PosEvidence(1, 1, point))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf pos" begin
    @test rf ≈ norm(point, Inf) - 1/ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ norm(point, Inf)
end

## Lie
nloc = 1
gen = CPB.Generator(nvar, nloc)
M = 8.0

CPB.add_af!(gen, 1)
CPB.add_evidence!(gen, NegEvidence(1, [0, 0]))
CPB.add_evidence!(gen, PosEvidence(1, 1, [8, 0]))
point1 = [2, 0]
point2 = [4, 0]
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 1, point2, 0.0))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf lie" begin
    @test rf ≈ (4 - δ - 1/ϵ)/2
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 4
end

## Lie
nloc = 1
gen = CPB.Generator(nvar, nloc)
M = 32.0

CPB.add_af!(gen, 1)
CPB.add_evidence!(gen, NegEvidence(1, [0, 0]))
CPB.add_evidence!(gen, PosEvidence(1, 1, [8, 0]))
point1 = [6, 0]
point2 = [4, 0]
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 1, point2, 0.0))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf lie" begin
    @test rf ≈ 5
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 6
end

## Pos and Lie: 2 wits #1
nloc = 2
gen = CPB.Generator(nvar, nloc)
M = 32.0

CPB.add_af!(gen, 1)
CPB.add_af!(gen, 2)
CPB.add_evidence!(gen, NegEvidence(1, [-2, 0]))
CPB.add_evidence!(gen, NegEvidence(1, [2, 0]))
CPB.add_evidence!(gen, PosEvidence(2, 1, [2, 0]))
point1 = [2, 0]
point2 = [4, 0]
nA = 0.5
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 2, point2, nA))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf pos and lie: 2 wits #1" begin
    @test rf ≈ (2 - δ - 1/ϵ)/2
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 2/(1 + nA)
end

## Pos and Lie: 2 wits #2
nloc = 2
gen = CPB.Generator(nvar, nloc)
M = 32.0

CPB.add_af!(gen, 1)
CPB.add_af!(gen, 2)
CPB.add_evidence!(gen, NegEvidence(1, [-2, 0]))
CPB.add_evidence!(gen, NegEvidence(1, [2, 0]))
CPB.add_evidence!(gen, PosEvidence(2, 1, [2, 0]))
point1 = [4, 0]
point2 = [2, 0]
nA = 0.5
CPB.add_evidence!(gen, LieEvidence(1, 1, point1, 2, point2, nA))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf pos and lie: 2 wits #2" begin
    @test rf ≈ 2 - δ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 2
end

## Pos and Lie: 2 wits #3
nloc = 2
gen = CPB.Generator(nvar, nloc)
M = 32.0

CPB.add_af!(gen, 1)
CPB.add_af!(gen, 1)
CPB.add_af!(gen, 2)
CPB.add_evidence!(gen, NegEvidence(1, [-2, 0]))
CPB.add_evidence!(gen, NegEvidence(1, [2, 0]))
CPB.add_evidence!(gen, PosEvidence(1, 1, [-4, 0]))
CPB.add_evidence!(gen, PosEvidence(1, 2, [4, 0]))
point1 = [1, 0]
point2 = [1, 0]
point3 = [2, 0]
nA = 1.0
CPB.add_evidence!(gen, LieEvidence(1, 2, point1, 2, point2, nA))
CPB.add_evidence!(gen, LieEvidence(2, 1, point2, 1, point3, nA))

ϵ = 1e5
δ = 1.0
rf = CPB.compute_mpf_feasibility(gen, ϵ, δ, M, solver)[2]
mpf, re = CPB.compute_mpf_evidence(gen, M, solver)

@testset "compute mpf pos and lie: 2 wits #3" begin
    @test rf ≈ (1 - 1/ϵ)/2
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
    @test re ≈ 1
end

nothing