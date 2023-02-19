module ExampleDai2021

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

mpf_inv = MultiPolyFunc([PolyFunc([
    AffForm([-1.0, 0.0], -2.0),
    AffForm([1.0, 0.0], -2.0),
    AffForm([0.0, -1.0], -2.0),
    AffForm([0.0, 1.0], -2.0)
])])

pf_dom = PolyFunc([
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, 1.0], 0.0)
])
A = [-0.999 0.0; -0.139 0.341]
b = [0.0, 0.0]
piece1 = Piece(pf_dom, 1, A, b, 1)
#
pf_dom = PolyFunc([
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
])
A = [0.436 0.323; 0.388 -0.049]
b = [0.0, 0.0]
piece2 = Piece(pf_dom, 1, A, b, 1)
#
pf_dom = PolyFunc([
    AffForm([1.0, 0.0], 0.0), AffForm([0.0, 1.0], 0.0)
])
A = [-0.457 0.215; 0.491 0.49]
b = [0.0, 0.0]
piece3 = Piece(pf_dom, 1, A, b, 1)
#
pf_dom = PolyFunc([
    AffForm([1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
])
A = [-0.022 0.344; 0.458 0.271]
b = [0.0, 0.0]
piece4 = Piece(pf_dom, 1, A, b, 1)
#
sys = System([piece1, piece2, piece3, piece4])

mpf_safe = MultiPolyFunc([PolyFunc([
    AffForm([-1.0, -1.0], -1.8),
    AffForm([-1.0, 1.0], -1.8),
    AffForm([1.0, -1.0], -1.8),
    AffForm([1.0, 1.0], -1.8)
])])

mlist_init = [[[-1.0, 0.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]]]

# Illustration
fig = figure(0, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=2, ncols=4,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("aspect"=>"equal")
)

xlims = (-3.2, 3.2)
ylims = (-3.2, 3.2)
lims = [(-10, -10), (10, 10)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.plot(0, 0, marker="x", ms=7, c="black", mew=1.5)

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
const wit_trace = Witness[]
const mpf_trace = MultiPolyFunc{PolyFunc{AffForm{Vector{Float64},Float64}}}[]
function rec_mpf_wit_trace(::Any, mpf, wit)
    push!(mpf_trace, MultiPolyFunc(
        map(pf -> PolyFunc(copy(pf.afs)), mpf.pfs
    )))
    push!(wit_trace, Witness(
        copy.(wit.mlist_inside), copy.(wit.mlist_image),
        copy.(wit.mlist_unknown), copy.(wit.mlist_outside)
    ))
end

ϵ = 0.1
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=false, callback_fcn=rec_mpf_wit_trace
)

display(status)

push!(mpf_trace, mpf)
push!(wit_trace, wit)

for (iter, (wit, mpf)) in enumerate(zip(wit_trace, mpf_trace))
    for (loc, points) in enumerate(wit.mlist_inside)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="blue", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.mlist_image)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="purple", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.mlist_outside)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="red", ms=5)
        end
    end

    for (loc, points) in enumerate(wit.mlist_unknown)
        @assert loc == 1
        for point in points
            plot_point!(ax_[iter], point, mc="orange", ms=5)
        end
    end

    for (loc, pf) in enumerate(mpf.pfs)
        @assert loc == 1
        plot_level!(ax_[iter], pf.afs, lims)
    end
end

fig.savefig(string(
    @__DIR__, "/../figures/fig_exa_dai2020_1.png"
), dpi=200, transparent=false, bbox_inches="tight")

end # module