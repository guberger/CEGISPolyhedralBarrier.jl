module ExampleSeeSaw

using LinearAlgebra
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
Piece = CPB.Piece
System = CPB.System
Witness = CPB.Witness

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

N = 2
M = 1

INN = Matrix{Float64}(I, N, N)

mpf_inv = MultiPolyFunc([PolyFunc([
    AffForm([-1.0, 0.0], -20.0),
    AffForm([1.0, 0.0], -20.0),
    AffForm([0.0, -1.0], -20.0),
    AffForm([0.0, 1.0], -20.0)
])])

pf_dom = PolyFunc([AffForm([1.0, 0.0], -4.0)])
b = [1.0, 2.0]
piece1 = Piece(pf_dom, 1, INN, b, 1)
#
pf_dom = PolyFunc([
    AffForm([-1.0, 0.0], 5.0), AffForm([1.0, 0.0], -7.0)
])
b = [2.0, 1.0]
piece2 = Piece(pf_dom, 1, INN, b, 1)
#
pf_dom = PolyFunc([
    AffForm([-1.0, 0.0], 7.0), AffForm([1.0, 0.0], -9.0)
])
b = [1.0, 3.0]
piece3 = Piece(pf_dom, 1, INN, b, 1)
#
pf_dom = PolyFunc([AffForm([-1.0, 0.0], 9.0)])
b = [2.0, 1.0]
piece4 = Piece(pf_dom, 1, INN, b, 1)
#
sys = System([piece1, piece2, piece3, piece4])

mlist_init = [[[0.0, 0.0]]]

mpf_safe = MultiPolyFunc([PolyFunc([AffForm([1.0, -2.1], -0.1)])])

# Illustration
fig = figure(0, figsize=(10, 8))
ax = fig.add_subplot(aspect="equal")
ax_ = (ax,)

xlims = (-22, 22)
ylims = (-22, 22)
lims = [(-40, -40), (40, 40)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.plot(0, 0, marker="x", ms=10, c="black", mew=2.5)
end

for (loc, pf) in enumerate(mpf_safe.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="green", fa=0.1, ec="green")
end

for (loc, pf) in enumerate(mpf_inv.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", ec="yellow")
end

for piece in sys.pieces
    plot_level!(
        ax_[piece.loc1], piece.pf_dom.afs, lims, fc="blue", fa=0.1, ec="blue"
    )
end

## Learner
ϵ = 1e-2
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver
)

display(status)

for (loc, pf) in enumerate(mpf.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="red", ec="red", fa=0.1, ew=0.5)
end

for (loc, points) in enumerate(wit.mlist_inside)
    for point in points
        plot_point!(ax_[loc], point, mc="blue")
    end
end

for (loc, points) in enumerate(wit.mlist_image)
    for point in points
        plot_point!(ax_[loc], point, mc="purple")
    end
end

for (loc, points) in enumerate(wit.mlist_outside)
    for point in points
        plot_point!(ax_[loc], point, mc="red")
    end
end

for (loc, points) in enumerate(wit.mlist_unknown)
    for point in points
        plot_point!(ax_[loc], point, mc="orange")
    end
end

fig.savefig(string(
    @__DIR__, "/../figures/fig_exa_see_saw.png"
), dpi=200, transparent=false, bbox_inches="tight")

end # module