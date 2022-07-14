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
PolyFunc = CPB.PolyFunc
_norm(pf::PolyFunc) =
    isempty(pf.afs) ? 0.0 : maximum(af -> norm(af.a, Inf), pf.afs)

solver() = Model(optimizer_with_attributes(
    HiGHS.Optimizer, "output_flag"=>false
))

ϵ = 1e-2
βmax = 100

## Empty
gen = CPB.Generator{2,1}()

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf empty" begin
    @test r ≈ 10
    @test all(pf -> isempty(pf.afs), mpf.pfs)
end

## Set #1
gen = CPB.Generator{2,1}()

CPB.add_neg_evidence!(gen, 1, SVector(0.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(0.5, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #1" begin
    @test r ≈ 0.5/2 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

## Set #2
gen = CPB.Generator{2,1}()

CPB.add_neg_evidence!(gen, 1, SVector(0.0, 0.0))
CPB.add_neg_evidence!(gen, 1, SVector(4.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(8.0, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #2" begin
    @test r ≈ 2 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

## Set #3
gen = CPB.Generator{2,1}()

CPB.add_neg_evidence!(gen, 1, SVector(0.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(6.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(8.0, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #3" begin
    @test r ≈ 3 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

## Set #4
gen = CPB.Generator{2,2}()

CPB.add_neg_evidence!(gen, 1, SVector(-2.0, 0.0))
CPB.add_neg_evidence!(gen, 1, SVector(2.0, 0.0))
CPB.add_neg_evidence!(gen, 2, SVector(4.0, 0.0))
CPB.add_pos_evidence!(gen, 2, SVector(2.0, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #4" begin
    @test r ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

## Set #5
gen = CPB.Generator{2,2}()

CPB.add_neg_evidence!(gen, 1, SVector(-2.0, 0.0))
CPB.add_neg_evidence!(gen, 1, SVector(2.0, 0.0))
CPB.add_pos_evidence!(gen, 2, SVector(2.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(4.0, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #5" begin
    @test r ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

## Set #6
gen = CPB.Generator{2,2}()

CPB.add_neg_evidence!(gen, 1, SVector(-2.0, 0.0))
CPB.add_neg_evidence!(gen, 1, SVector(2.0, 0.0))
CPB.add_neg_evidence!(gen, 2, SVector(1.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(-4.0, 0.0))
CPB.add_pos_evidence!(gen, 1, SVector(4.0, 0.0))

mpf, r = CPB.compute_mpf_robust(gen, ϵ, βmax, solver)

@testset "compute mpf #6" begin
    @test r ≈ 1 - ϵ
    @test maximum(pf -> _norm(pf), mpf.pfs) ≈ 1
end

nothing