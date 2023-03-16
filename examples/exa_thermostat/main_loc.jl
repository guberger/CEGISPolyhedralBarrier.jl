module ExampleThermostat

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
const _AT = AffForm{Vector{Float64},Float64}
const _PT = PolyFunc{_AT}

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

N = 2
M = 4

Tlo = 0.0
Ilo = 15.0
Slo = 18.0
Llo = 19.0
Lup = 21.0
Sup = 22.0
Iup = 25.0
Tup = 30.0
cf = 1.0
dt = 0.03
Tstab = 1.5 - dt

locs = [[1, 3], [2, 4]]
guards_temp = [
    [AffForm([-1.0, 0.0], Llo), AffForm([1.0, 0.0], -Llo)],
    [AffForm([-1.0, 0.0], Lup), AffForm([1.0, 0.0], -Lup)]
]
guards_time = [AffForm([0.0, 1.0], -Tstab), AffForm([0.0, -1.0], Tstab)]
α = exp(-cf*dt)
As = [[α 0; 0 1], [α 0; 0 0]]
bs = [
    [[Tlo*(1 - α), dt], [Tlo*(1 - α), Tstab + dt]],
    [[Tup*(1 - α), dt], [Tup*(1 - α), Tstab + dt]]
]

pieces = Piece{_PT,Matrix{Float64},Vector{Float64}}[]
for (i1, j1, i2, j2) in Iterators.product(1:2, 1:2, 1:2, 1:2)
    j1 > j2 && continue
    pf_dom = PolyFunc([guards_temp[i1][i2], guards_time[j2]])
    push!(pieces, Piece(pf_dom, locs[i1][j1], As[j2], bs[i2][j2], locs[i2][j2]))
end
sys = System(pieces)

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

mpf_inv = MultiPolyFunc([PolyFunc(_AT[]) for loc = 1:4])

mpf_safe = MultiPolyFunc([
    [PolyFunc([
        AffForm([-1.0, 0.0], Tlo), AffForm([1.0, 0.0], -Tup),
        AffForm([0.0, -1.0], 0.0), AffForm([0.0, 1.0], -(Tstab + 2*dt))
    ]) for loc = 1:2]...,
    [PolyFunc([
        AffForm([-1.0, 0.0], Slo), AffForm([1.0, 0.0], -Sup),
        AffForm([0.0, -1.0], 0.0), AffForm([0.0, 1.0], -(Tstab + 2*dt))
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

ϵ = dt/3
display(ϵ)
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=true
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

f = open(string(@__DIR__, "/witnesses.txt"), "w")
println(f, "image")
for (loc, points_image) in enumerate(wit.mlist_image)
    str = string(loc, ":")
    for point in points_image
        str = string(str, "(", point[1], ",", point[2], ")")
    end
    println(f, str)
end
println(f, "unknown")
for (loc, points_unknown) in enumerate(wit.mlist_unknown)
    str = string(loc, ":")
    for point in points_unknown
        str = string(str, "(", point[1], ",", point[2], ")")
    end
    println(f, str)
end
println(f, "outside")
for (loc, points_outside) in enumerate(wit.mlist_outside)
    str = string(loc, ":")
    for point in points_outside
        str = string(str, "(", point[1], ",", point[2], ")")
    end
    println(f, str)
end
println(f, "invariant")
zip_iter = zip(mpf_safe.pfs, mpf_inv.pfs, mpf.pfs)
for (loc, (pf_safe, pf_inv, pf)) in enumerate(zip_iter)
    str = string(loc, ":")
    afs = union(pf_safe.afs, pf_inv.afs, pf.afs)
    A = zeros(length(afs), 2)
    b = zeros(length(afs))
    for (i, af) in enumerate(afs)
        A[i, 1], A[i, 2] = (af.a...,)
        b[i] = -af.β
    end
    verts = compute_vertices_hrep(A, b)
    for vert in verts
        str = string(str, "(", vert[1], ",", vert[2], ")")
    end
    if !isempty(verts)
        str = string(str, "(", verts[1][1], ",", verts[1][2], ")")
    end
    println(f, str)
end
close(f)

end # module