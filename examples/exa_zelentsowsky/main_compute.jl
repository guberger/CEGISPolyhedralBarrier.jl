module ExampleZelentsowsky_Compute

using LinearAlgebra
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier

const GUROBI_ENV = Gurobi.Env()
function solver()
    model = direct_model(
        optimizer_with_attributes(() -> Gurobi.Optimizer(GUROBI_ENV))
    )
    set_optimizer_attribute(model, "OutputFlag", 0)
    # set_optimizer_attribute(model, "Method", 2)
    return model
end

datafile = "dataset_1"
include(string("./datasets/", datafile, ".jl"))

## Parameters
nvar = 2
nloc = 1

sys = CPB.System()

domain = CPB.Cone()
A = [0.0 1.0; -2.0 -1.0]
CPB.add_piece_cont!(sys, domain, 1, A)

domain = CPB.Cone()
B = [0.0 0.0; -1.0 0.0]
CPB.add_piece_cont!(sys, domain, 1, A + α*B)

## Learner feasible illustration
lear = CPB.Learner(nvar, nloc, sys, τ, ϵ, δ)
CPB.set_tol!(lear, :rad, 1e-6)

points_init = [[-1.0, 0.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]]
for point in points_init
    CPB.add_witness!(lear, 1, point)
end

## Solving
status, mpf, niter = CPB.learn_lyapunov!(lear, 1000, solver, solver)

display(status)

f = open(string(@__DIR__, "/results/", datafile, ".txt"), "w")
for lf in mpf.pfs[1].lfs
    println(f, lf.lin)
end
close(f)

end # module