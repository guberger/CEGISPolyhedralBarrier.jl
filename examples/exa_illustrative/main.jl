module ExampleIllustrative

using LinearAlgebra
using StaticArrays
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
System = CPB.System
MultiSet = CPB.MultiSet
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

mpf_inv = MultiPolyFunc{2,2}()
for loc = 1:2
    CPB.add_af!(mpf_inv, loc, SVector(-1.0, 0.0), -2.0)
    CPB.add_af!(mpf_inv, loc, SVector(1.0, 0.0), -2.0)
    CPB.add_af!(mpf_inv, loc, SVector(0.0, -1.0), -2.0)
    CPB.add_af!(mpf_inv, loc, SVector(0.0, 1.0), -2.0)
end

sys = System{2}()

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(0.0, -1.0), 0.5)
A = @SMatrix [0.5 0.0; 0.0 0.5]
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 2)

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(1.0, 0.0), 0.0)
A = @SMatrix [1.0 0.0; 0.0 1.0]
b = @SVector [0.0, 0.5]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

mset_init = MultiSet{2,2}()
CPB.add_point!(mset_init, 1, SVector(-1.0, -1.0))
CPB.add_point!(mset_init, 1, SVector(-1.0, 1.0))
CPB.add_point!(mset_init, 1, SVector(1.0, -1.0))
CPB.add_point!(mset_init, 1, SVector(1.0, 1.0))

mpf_safe = MultiPolyFunc{2,2}()
CPB.add_af!(mpf_safe, 1, SVector(0.0, 1.0), -1.5)

# Illustration
fig = figure(0, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=1, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("aspect"=>"equal")
)

xlims = (-2.2, 2.2)
ylims = (-2.2, 2.2)
lims = [(-10, -10), (10, 10)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.tick_params(axis="both", labelsize=15)
    ax.plot(0, 0, marker="x", ms=10, c="black", mew=2.5)
end

for (loc, pf) in enumerate(mpf_safe.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", fa=0, ec="red", ew=1.5)
end

for (loc, pf) in enumerate(mpf_inv.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", ec="blue")
end

# for piece in sys.pieces
#     plot_level!(
#         ax_[piece.loc1], piece.pf_dom.afs, lims, fc="blue", fa=0.1, ec="blue"
#     )
# end

## Learner
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, 1e-1, 1e-8)
status, mpf, wit = CPB.learn_lyapunov!(lear, Inf, solver, solver)

display(status)

for (loc, pf) in enumerate(mpf.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="gold", ec="gold", fa=0.5, ew=2.5)
end

for (loc, points) in enumerate(wit.inside.sets)
    for point in points
        plot_point!(ax_[loc], point, mc="blue")
    end
end

for (loc, points) in enumerate(wit.image.sets)
    for point in points
        plot_point!(ax_[loc], point, mc="purple")
    end
end

for (loc, points) in enumerate(wit.outside.sets)
    for point in points
        plot_point!(ax_[loc], point, mc="red")
    end
end

for (loc, points) in enumerate(wit.unknown.sets)
    for point in points
        plot_point!(ax_[loc], point, mc="orange")
    end
end

fig.savefig(string(
    @__DIR__, "/../figures/fig_exa_illustrative.png"
), dpi=200, transparent=false, bbox_inches="tight")

end # module