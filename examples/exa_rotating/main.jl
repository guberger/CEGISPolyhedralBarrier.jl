module ExampleRotating

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

mpf_safe = MultiPolyFunc{2,2}()
for loc = 1:2
    CPB.add_af!(mpf_safe, loc, SVector(-1.0, 0.0), -2.0)
    CPB.add_af!(mpf_safe, loc, SVector(1.0, 0.0), -2.0)
    CPB.add_af!(mpf_safe, loc, SVector(0.0, -1.0), -2.0)
    CPB.add_af!(mpf_safe, loc, SVector(0.0, 1.0), -2.0)
end

sys = System{2}()
θ = π/5
α = 0.9

pf_dom = PolyFunc{2}()
A = α*(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 1, A, b, 2)

pf_dom = PolyFunc{2}()
A = α*(@SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)])
b = @SVector [0.0, 0.0]
CPB.add_piece!(sys, pf_dom, 2, A, b, 1)

mset_init = MultiSet{2,2}()
CPB.add_point!(mset_init, 1, SVector(-0.25, -0.25))
CPB.add_point!(mset_init, 1, SVector(-0.25, 0.25))
CPB.add_point!(mset_init, 1, SVector(0.25, -0.25))
CPB.add_point!(mset_init, 1, SVector(0.25, 0.25))

mpf_inv = MultiPolyFunc{2,2}()

## Learner
lear = CPB.Learner(sys, mpf_safe, mpf_inv, mset_init, 0.1, 1e-8)
const _mpfs = MultiPolyFunc{2,2}[]
const _wits = CPB.Witness{2,2}[]
function rec_mpfs_wits(::Any, mpf, wit)
    push!(_mpfs, CPB.MultiPolyFunc(
        map(pf -> PolyFunc(copy(pf.afs)), mpf.pfs)
    ))
    inside_ = MultiSet(copy.(wit.inside.sets))
    image_ = MultiSet(copy.(wit.image.sets))
    outside_ = MultiSet(copy.(wit.outside.sets))
    unknown_ = MultiSet(copy.(wit.unknown.sets))
    push!(_wits, CPB.Witness(inside_, image_, outside_, unknown_))
end
status, mpf, wit = CPB.learn_lyapunov!(
    lear, Inf, solver, solver, callback_fcn=rec_mpfs_wits
)

display(status)

# Illustration
fig = figure(0, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=1, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("aspect"=>"equal")
)

xlims = (-2.8, 2.8)
ylims = (-2.8, 2.8)
lims = [(-10, -10), (10, 10)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.tick_params(axis="both", labelsize=15)
    ax.plot(0, 0, marker="x", ms=10, c="black", mew=2.5)
end

for (loc, pf) in enumerate(mpf_safe.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", fa=0, ec="green")
end

hinsides = map(ax -> ax.plot((), ())[1], ax_)
for h in hinsides
    h.set_linestyle("none")
    h.set_color("blue")
    h.set_marker(".")
    h.set_markersize(10)
end
himages = map(ax -> ax.plot((), ())[1], ax_)
for h in himages
    h.set_linestyle("none")
    h.set_color("purple")
    h.set_marker(".")
    h.set_markersize(10)
end
houtsides = map(ax -> ax.plot((), ())[1], ax_)
for h in houtsides
    h.set_linestyle("none")
    h.set_color("red")
    h.set_marker(".")
    h.set_markersize(10)
end
hunknowns = map(ax -> ax.plot((), ())[1], ax_)
for h in hunknowns
    h.set_linestyle("none")
    h.set_color("orange")
    h.set_marker(".")
    h.set_markersize(10)
end
polypfs = map(ax -> matplotlib.collections.PolyCollection((zeros(0, 2),)), ax_)
for (ax, p) in zip(ax_, polypfs)
    p.set_facecolor("none")
    p.set_edgecolor("red")
    p.set_linewidth(1.0)
    ax.add_collection(p)
end

for (iter, (wit, mpf)) in enumerate(zip(_wits, _mpfs))
    for (h, points) in zip(hinsides, wit.inside.sets)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(himages, wit.image.sets)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(houtsides, wit.outside.sets)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(hunknowns, wit.unknown.sets)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (p, pf) in zip(polypfs, mpf.pfs)
        A = zeros(length(pf.afs), 2)
        b = zeros(length(pf.afs))
        for (i, af) in enumerate(pf.afs)
            A[i, 1], A[i, 2] = (af.a...,)
            b[i] = -af.β
        end
        A = vcat(A, -Matrix{Bool}(I, 2, 2), Matrix{Bool}(I, 2, 2))
        b = vcat(b, -collect(lims[1]), collect(lims[2]))
        verts = compute_vertices_hrep(A, b)
        # isempty(verts) && continue
        p.set_verts((verts,))
    end
    fig.canvas.draw()
    fig.canvas.flush_events()
    fig.savefig(string(
        "./examples/figures/animation_rotating/frame_", iter, ".png"
    ), bbox_inches="tight", dpi=50)
end

end # module