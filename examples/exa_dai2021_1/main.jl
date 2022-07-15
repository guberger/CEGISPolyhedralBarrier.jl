module ExampleDai2021

using LinearAlgebra
using StaticArrays
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
System = CPB.System
PointSet = CPB.PointSet
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

mpf_inv = MultiPolyFunc{2,1}()
CPB.add_af!(mpf_inv, 1, SVector(-1.0, 0.0), -2.0)
CPB.add_af!(mpf_inv, 1, SVector(1.0, 0.0), -2.0)
CPB.add_af!(mpf_inv, 1, SVector(0.0, -1.0), -2.0)
CPB.add_af!(mpf_inv, 1, SVector(0.0, 1.0), -2.0)

sys = System{2}()

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(-1.0, 0.0), 0.0)
CPB.add_af!(pf_dom, SVector(0.0, 1.0), 0.0)
A = @SMatrix [-0.999 0.0; -0.139 0.341]
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(-1.0, 0.0), 0.0)
CPB.add_af!(pf_dom, SVector(0.0, -1.0), 0.0)
A = @SMatrix [0.436 0.323; 0.388 -0.049]
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(1.0, 0.0), 0.0)
CPB.add_af!(pf_dom, SVector(0.0, 1.0), 0.0)
A = @SMatrix [-0.457 0.215; 0.491 0.49]
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

pf_dom = PolyFunc{2}()
CPB.add_af!(pf_dom, SVector(1.0, 0.0), 0.0)
CPB.add_af!(pf_dom, SVector(0.0, -1.0), 0.0)
A = @SMatrix [-0.022 0.344; 0.458 0.271]
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 1)

mpf_safe = MultiPolyFunc{2,1}()
CPB.add_af!(mpf_safe, 1, SVector(-1.0, -1.0), -1.8)
CPB.add_af!(mpf_safe, 1, SVector(-1.0, 1.0), -1.8)
CPB.add_af!(mpf_safe, 1, SVector(1.0, -1.0), -1.8)
CPB.add_af!(mpf_safe, 1, SVector(1.0, 1.0), -1.8)

iset = PointSet{2,1}()
CPB.add_point!(iset, 1, SVector(-1.0, 0.0))
CPB.add_point!(iset, 1, SVector(1.0, 0.0))
CPB.add_point!(iset, 1, SVector(0.0, -1.0))
CPB.add_point!(iset, 1, SVector(0.0, 1.0))

# Illustration
fig = figure(0, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=3, ncols=4,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("aspect"=>"equal")
)

xlims = (-4.2, 4.2)
ylims = (-4.2, 4.2)
lims = [(-10, -10), (10, 10)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.plot(0, 0, marker="x", ms=10, c="black", mew=2.5)

    for (loc, pf) in enumerate(mpf_safe.pfs)
        @assert loc == 1
        plot_level!(ax, pf.afs, lims, fc="green", fa=0.1, ec="green")
    end

    for (loc, pf) in enumerate(mpf_inv.pfs)
        @assert loc == 1
        plot_level!(ax, pf.afs, lims, fc="none", ec="yellow")
    end

    for piece in sys.pieces
        @assert piece.loc1 == 1
        plot_level!(ax, piece.pf_dom.afs, lims, fc="blue", fa=0.1, ec="blue")
    end
end


## Learner
lear = CPB.Learner(sys, mpf_safe, mpf_inv, iset, 0.1, 1e-8)
rec = CPB.TraceRecorder(lear)
status, mpf, wit = CPB.learn_lyapunov!(lear, Inf, solver, solver, rec=rec)

display(status)

for (loc, pf) in enumerate(mpf.pfs)
    plot_level!(
        ax_[length(rec.wit_list) + 1], pf.afs, lims,
        fc="red", ec="red", fa=0.1, ew=0.5
    )
end

for (iter, wit) in enumerate(rec.wit_list)
    for (loc, points) in enumerate(wit.soft_evid.points_list)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="purple", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.hard_evid.points_list)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="blue", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.unsafe.points_list)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="red", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.pos.points_list)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="orange", ms=5)
        end
    end
end

fig.savefig(string(
    @__DIR__, "/../figures/fig_exa_dai2020_1.png"
), dpi=200, transparent=false, bbox_inches="tight")

end # module