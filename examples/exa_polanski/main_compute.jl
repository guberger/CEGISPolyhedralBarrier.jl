module ExamplePolanski_Compute

using LinearAlgebra
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier

const GUROBI_ENV = Gurobi.Env()
solver = optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
)

datafile = "dataset_2"
include(string("./datasets/", datafile, ".jl"))

## Parameters
nvar = 3

sys = CPB.System()

domain = CPB.Cone()
A = [-10.0 -2.0 -2.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
CPB.add_piece!(sys, domain, A)

domain = CPB.Cone()
A = [-10.0 -10.0*α -10.0*α; 1.0 0.0 0.0; 0.0 1.0 0.0]
CPB.add_piece!(sys, domain, A)

## Learner feasible illustration
lear = CPB.Learner(nvar, sys, ϵ, θ, δ)
CPB.set_tol!(lear, :rad, 1e-6)

points_init = [
    [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0],
    [0.0, -1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, -1.0], [0.0, 0.0, 1.0]
]
for point in points_init
    CPB.add_point_init!(lear, point)
end

## Solving
sol = CPB.learn_lyapunov!(lear, 1000, solver)

display(sol.status)

# f = open(string(@__DIR__, "/results/", datafile, ".txt"), "w")
# for lf in sol.lfs_list[sol.niter]
#     println(f, lf.lin)
# end
# close(f)

end # module