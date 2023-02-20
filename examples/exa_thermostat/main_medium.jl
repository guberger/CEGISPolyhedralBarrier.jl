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
M = 4

Tlo = 0
Ilo = 15
Slo = 18
Llo = 19
Lup = 21
Sup = 22
Iup = 25
Tup = 30
cf = 1.0
dt = 0.03
Tstab = 1.5

A = [exp(-cf*dt) 0; 0 1]
Astab = [exp(-cf*dt) 0; 0 0]
blo = [Tlo*(1 - exp(-cf*dt)), dt]
blostab = [Tlo*(1 - exp(-cf*dt)), 0]
bup = [Tup*(1 - exp(-cf*dt)), dt]
bupstab = [Tup*(1 - exp(-cf*dt)), 0]

# cooling
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Llo), AffForm([0.0, 1.0], -Tstab)])
piece1a = Piece(pf_dom, 1, A, blo, 1)
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Llo), AffForm([0.0, -1.0], Tstab)])
piece1b = Piece(pf_dom, 1, Astab, blostab, 3)
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Llo)])
piece1c = Piece(pf_dom, 3, Astab, blostab, 3)
# start heating
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Llo), AffForm([0.0, 1.0], -Tstab)])
piece2a = Piece(pf_dom, 1, A, bup, 2)
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Llo), AffForm([0.0, -1.0], Tstab)])
piece2b = Piece(pf_dom, 1, Astab, bupstab, 4)
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Llo)])
piece2c = Piece(pf_dom, 3, Astab, bupstab, 4)
# heating
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Lup), AffForm([0.0, 1.0], -Tstab)])
piece3a = Piece(pf_dom, 2, A, bup, 2)
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Lup), AffForm([0.0, -1.0], Tstab)])
piece3b = Piece(pf_dom, 2, Astab, bupstab, 4)
pf_dom = PolyFunc([AffForm([1.0, 0.0], -Lup)])
piece3c = Piece(pf_dom, 4, Astab, bupstab, 4)
# start cooling
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Lup), AffForm([0.0, 1.0], -Tstab)])
piece4a = Piece(pf_dom, 2, A, blo, 1)
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Lup), AffForm([0.0, -1.0], Tstab)])
piece4b = Piece(pf_dom, 2, Astab, blostab, 3)
pf_dom = PolyFunc([AffForm([-1.0, 0.0], Lup)])
piece4c = Piece(pf_dom, 4, Astab, blostab, 3)
#
sys = System([
    piece1a, piece1b, piece1c,
    piece2a, piece2b, piece2c,
    piece3a, piece3b, piece3c,
    piece4a, piece4b, piece4c
])

# simulation
fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 300
ax.set_xlim((0, nstep*dt))
ax.set_xlabel("time")

DT = Tup - Tlo
ax.set_ylim((Tlo - 0.05*DT, Tup + 0.05*DT))
for y in (Tlo, Slo, Llo, Lup, Sup, Tup)
    ax.plot((0, nstep), (y, y), c="k")
end
ax.set_ylabel("temperature")

T0 = Tlo
x = [T0, 0.0]
loc = 2
tol_dom = 1e-8

for istep = 1:nstep
    global x, loc
    @assert loc ∈ (1, 2, 3, 4)
    color = loc ∈ (1, 3) ? "blue" : "red"
    marker = loc ∈ (1, 2) ? "x" : "."
    ax.plot((istep - 1)*dt, x[1], marker=marker, ms=5, c=color)
    istep == nstep && break
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !CPB._prox(piece.pf_dom, x, 0) && continue
        loc = piece.loc2
        x = piece.A*x + piece.b
        break
    end
end

mpf_inv = MultiPolyFunc([
    PolyFunc(AffForm{Vector{Float64},Float64}[]) for loc = 1:4
])

mpf_safe = MultiPolyFunc([
    [PolyFunc([
        AffForm([-1.0, 0.0], Tlo), AffForm([1.0, 0.0], -Tup),
        AffForm([0.0, -1.0], 0), AffForm([0.0, 1.0], -(Tstab + 3*dt))
    ]) for loc = 1:2]...,
    [PolyFunc([
        AffForm([-1.0, 0.0], Slo), AffForm([1.0, 0.0], -Sup),
        AffForm([0.0, -1.0], -2*dt), AffForm([0.0, 1.0], -2*dt)
    ]) for loc = 1:2]...
])
# empty!(mpf_safe.pfs[3].afs)
# empty!(mpf_safe.pfs[4].afs)

mlist_init = [
    [[Ilo, 0.0], [Iup, 0.0]],
    [[Ilo, 0.0], [Iup, 0.0]],
    Vector{Float64}[],
    Vector{Float64}[]
]

ϵ = 0.005
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=false
)

display(status)

# Illustration
fig = figure(1, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=2, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1)
)

xlims = (Tlo - 0.05*DT, Tup + 0.05*DT)
ylims = (-5*dt, Tstab + 5*dt)
lims = ([Tlo - 2*DT, -Tstab - 5*dt], [Tup + 2*DT, 2*Tstab + 5*dt])

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.tick_params(axis="both", labelsize=15)
end

for (loc, pf) in enumerate(mpf_safe.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", fa=0, ec="red", ew=1.5)
end

for (loc, pf) in enumerate(mpf_inv.pfs)
    plot_level!(ax_[loc], pf.afs, lims, fc="none", ec="blue")
end

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

end # module