module ExampleIllustrative

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
M = 2

mpf_inv = MultiPolyFunc([PolyFunc([
    AffForm([-1.0, 0.0], -2.0),
    AffForm([1.0, 0.0], -2.0),
    AffForm([0.0, -1.0], -2.0),
    AffForm([0.0, 1.0], -2.0)
]) for loc = 1:2])

pf_dom = PolyFunc([AffForm([0.0, -1.0], 0.5)])
A = [0.5 0.0; 0.0 0.5]
b = [0.0, 0.0]
piece1 = Piece(pf_dom, 1, A, b, 2)
pf_dom = PolyFunc([AffForm([1.0, 0.0], 0.0)])
A = [1.0 0.0; 0.0 1.0]
b = [0.0, 0.5]
piece2 = Piece(pf_dom, 2, A, b, 1)
sys = System([piece1, piece2])

mlist_init = [
    [[-1.0, -1.0], [-1.0, 1.0], [1.0, -1.0], [1.0, 1.0]], Vector{Float64}[]
]

mpf_safe = MultiPolyFunc([
    PolyFunc([AffForm([0.0, 1.0], -1.5)]),
    PolyFunc(AffForm{Vector{Float64},Float64}[])
])

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
ϵ = 0.1
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver
)

display(status)

for (loc, pf) in enumerate(mpf.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="gold", ec="gold", fa=0.5, ew=2.5)
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
    @__DIR__, "/../figures/fig_exa_illustrative.png"
), dpi=200, transparent=false, bbox_inches="tight")

end # module