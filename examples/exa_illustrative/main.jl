module ExampleIllustrative

using LinearAlgebra
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
Polyhedron = CPB.Polyhedron
PolyFunc = CPB.PolyFunc
System = CPB.System
InitSet = CPB.InitSet
UnsafeSet = CPB.UnsafeSet
State = CPB.State
Region = CPB.Region

include("../utils/geometry.jl")
include("plotting.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

## Parameters
nvar = 2
nloc = 2

box = Polyhedron()
CPB.add_halfspace!(box, [-1, 0], -2)
CPB.add_halfspace!(box, [1, 0], -2)
CPB.add_halfspace!(box, [0, -1], -2)
CPB.add_halfspace!(box, [0, 1], -2)

sys = System()

domain = Polyhedron()
CPB.add_halfspace!(domain, [0, -1], 0.5)
A = [0.5 0.0; 0.0 0.5]
b = [0.0, 0.0]
CPB.add_piece!(sys, domain ∩ box, 1, A, b, 2)

domain = Polyhedron()
CPB.add_halfspace!(domain, [1, 0], 0)
A = Matrix{Bool}(I, 2, 2)
b = [0.0, 0.5]
CPB.add_piece!(sys, domain ∩ box, 2, A, b, 1)

iset = InitSet()
init_points = ([-1, -1], [-1, 1], [1, -1], [1, 1])
for point in init_points
    CPB.add_state!(iset, 1, point)
end

uset = UnsafeSet()
udom = Polyhedron()
CPB.add_halfspace!(udom, [-1, 0], -2)
CPB.add_halfspace!(udom, [1, 0], 1)
CPB.add_halfspace!(udom, [0, -1], 1)
CPB.add_halfspace!(udom, [0, 1], -2)
CPB.add_region!(uset, 2, udom)

## Plotting

# Illustration
fig = figure(0, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=1, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("aspect"=>"equal")
)

xlims = (-2.2, 2.2)
ylims = (-2.2, 2.2)

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.set_xticks(-2:1:2)
    ax.set_yticks(-2:1:2)
    ax.tick_params(axis="both", labelsize=15)
    ax.plot(0, 0, marker="x", ms=10, c="black", mew=2.5)
end

idom = Polyhedron()
CPB.add_halfspace!(idom, [-1, 0], -1)
CPB.add_halfspace!(idom, [1, 0], -1)
CPB.add_halfspace!(idom, [0, -1], -1)
CPB.add_halfspace!(idom, [0, 1], -1)

for state in iset.states
    plot_point!(ax_[state.loc], state.point, mc="gold")
end
for loc = 1:nloc
    points = [state.point for state in iset.states if state.loc == loc]
    plot_vrep!(ax_[loc], points, fc="yellow", ec="yellow")
end

for region in uset.regions
    A, b = _to_vector_ineq(region.domain, 2)
    plot_hrep!(ax_[region.loc], A, b, fc="red", ec="red")
end

for piece in sys.pieces
    A, b = _to_vector_ineq(piece.domain, 2)
    plot_hrep!(ax_[piece.loc1], A, b, fa=0.1, ec="none")
end

## Learner
lear = CPB.Learner(nvar, nloc, sys, iset, uset, 0, 0)
CPB.set_tol!(lear, :rad, 1e-4)
CPB.set_tol!(lear, :bigM, 1e3)

status, mpf = CPB.learn_lyapunov!(lear, 1000, solver, solver)

display(status)

for loc = 1:nloc
    Abox, bbox = _to_vector_ineq(box, 2)
    Apf, bpf = _to_vector_ineq(mpf.pfs[loc], 2)
    plot_hrep!(ax_[loc], vcat(Abox, Apf), vcat(bbox, bpf))
end

end # module