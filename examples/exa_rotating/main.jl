module ExampleRotating

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

mpf_safe = MultiPolyFunc([PolyFunc([
    AffForm([-1.0, 0.0], -2.0),
    AffForm([1.0, 0.0], -2.0),
    AffForm([0.0, -1.0], -2.0),
    AffForm([0.0, 1.0], -2.0)
]) for loc = 1:2])

θ = π/5
α = 0.9

pf_dom = PolyFunc(AffForm{Vector{Float64},Float64}[])
A = α*[cos(θ) -sin(θ); sin(θ) cos(θ)]
b = [0.0, 0.0]
piece1 = Piece(pf_dom, 1, A, b, 2)
pf_dom = PolyFunc(AffForm{Vector{Float64},Float64}[])
A = α*[cos(θ) -sin(θ); sin(θ) cos(θ)]
b = [0.0, 0.0]
piece2 = Piece(pf_dom, 2, A, b, 1)
sys = System([piece1, piece2])

ρ = 0.25
mlist_init = [
    [[-ρ, -ρ], [-ρ, ρ], [ρ, -ρ], [ρ, ρ]], Vector{Float64}[]
]

mpf_inv = MultiPolyFunc([
    PolyFunc(AffForm{Vector{Float64},Float64}[]) for loc = 1:2
])

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

ϵ = 0.15
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=false, callback_fcn=rec_mpf_wit_trace
)

display(status)

push!(mpf_trace, mpf)
push!(wit_trace, wit)

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

niter = length(wit_trace)
nchar = 60
println(repeat("=", nchar))
ichar = 1

for (iter, (wit, mpf)) in enumerate(zip(wit_trace, mpf_trace))
    global ichar
    while iter ≥ ichar*niter/nchar
        ichar += 1
        print("=")
    end

    for (h, points) in zip(hinsides, wit.mlist_inside)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(himages, wit.mlist_image)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(houtsides, wit.mlist_outside)
        h.set_xdata(getindex.(points, 1))
        h.set_ydata(getindex.(points, 2))
    end
    for (h, points) in zip(hunknowns, wit.mlist_unknown)
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

println()

end # module