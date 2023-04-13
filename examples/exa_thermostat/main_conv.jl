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

N = 6
M = 1

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

γ = 100.0
hots = [[ifelse(j == i, 1.0, 0.0) for j = 1:6] for i = 1:6]
locs = [
    [[ifelse(k == 2*j - 2 + i, γ, 0.0) for k = 1:4] for j = 1:2] for i = 1:2
]
guards_temp = [
    [AffForm(-hots[1], +Llo), AffForm(+hots[1], -Llo)],
    [AffForm(-hots[1], +Lup), AffForm(+hots[1], -Lup)]
]
guards_time = [AffForm(+hots[2], -Tstab), AffForm(-hots[2], +Tstab)]
α = exp(-cf*dt)
As_ = [[α 0; 0 1], [α 0; 0 0]]
As = [[As_[i] zeros(2, 4); zeros(4, 6)] for i = 1:2]
bs_ = [
    [[Tlo*(1 - α), dt], [Tlo*(1 - α), Tstab + dt]],
    [[Tup*(1 - α), dt], [Tup*(1 - α), Tstab + dt]]
]
bs = [[[bs_[i][j]..., locs[i][j]...] for j = 1:2] for i = 1:2]
pieces = Piece{_PT,Matrix{Float64},Vector{Float64}}[]
for (i1, j1, i2, j2) in Iterators.product(1:2, 1:2, 1:2, 1:2)
    j1 > j2 && continue
    afs = [guards_temp[i1][i2], guards_time[j2]]
    loc1 = locs[i1][j1]
    for i = 3:6
        push!(afs, AffForm(+hots[i], -loc1[i - 2]))
        push!(afs, AffForm(-hots[i], +loc1[i - 2]))
    end
    pf_dom = PolyFunc(afs)
    push!(pieces, Piece(pf_dom, 1, As[j2], bs[i2][j2], 1))
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
x = [T0, 0.0, locs[1][1]...]
tol_dom = 1e-8

for istep = 1:nstep
    global x
    T = x[1]
    loc = x[3:6]
    @assert any(i -> any(j -> norm(loc - locs[i][j]) < 1e-5, 1:2), 1:2)
    color = any(j -> norm(loc - locs[1][j]) < 1e-5, 1:2) ? "blue" : "red"
    marker = any(i -> norm(loc - locs[i][1]) < 1e-5, 1:2) ? "x" : "."
    ax.plot((istep - 1)*dt, T, marker=marker, ms=5, c=color)
    istep == nstep && break
    flag = false
    for piece in sys.pieces
        @assert piece.loc1 == 1
        !CPB._prox(piece.pf_dom, x, 1e-9) && continue
        @assert piece.loc2 == 1
        flag = true
        x = piece.A*x + piece.b
        break
    end
    @assert flag
end

mpf_inv = MultiPolyFunc([PolyFunc(AffForm{Vector{Float64},Float64}[])])

mpf_safe = MultiPolyFunc([PolyFunc([
    AffForm([-1, 0, 0, 0, (Slo - Tlo)/γ, (Slo - Tlo)/γ], +Tlo),
    AffForm([+1, 0, 0, 0, (Tup - Sup)/γ, (Tup - Sup)/γ], -Tup),
    AffForm([0, -1, 0, 0, 0, 0], +0),
    AffForm([0, +1, 0, 0, 0, 0], -(Tstab + 2*dt))
])])
# empty!(mpf_safe.pfs[3].afs)
# empty!(mpf_safe.pfs[4].afs)

inits = [[Ilo, 0.0], [Iup, 0.0]]
mlist_init = [Vector{Float64}[]]
for (i, k) in Iterators.product(1:2, 1:2)
    push!(mlist_init[1], [inits[k]..., locs[i][1]...])
end

ϵ = dt/3
display(ϵ)
δ = 1e-8
iter_max = Inf

status, mpf, wit = @time CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=false, βmax=0.0
)

display(wit.mlist_inside[1])
display(wit.mlist_image[1])
display(wit.mlist_unknown[1])
display(wit.mlist_outside[1])
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

afs = AffForm{Vector{Float64},Float64}[]
locs_ = [locs[1]..., locs[2]...]

for (iloc, loc) in enumerate(locs_)
    empty!(afs)
    for af in mpf_safe.pfs[1].afs
        push!(afs, AffForm(af.a[1:2]*1.0, af.β + dot(af.a[3:6], loc*1.0)))
    end
    plot_level!(ax_[iloc], afs, lims, fc="none", fa=0, ec="red", ew=1.5)
    empty!(afs)
    for af in mpf_inv.pfs[1].afs
        push!(afs, AffForm(af.a[1:2]*1.0, af.β + dot(af.a[3:6], loc*1.0)))
    end
    plot_level!(ax_[iloc], afs, lims, fc="none", ec="blue")
    empty!(afs)
    for af in mpf.pfs[1].afs
        push!(afs, AffForm(af.a[1:2]*1.0, af.β + dot(af.a[3:6], loc*1.0)))
    end
    plot_level!(ax_[iloc], afs, lims, fc="gold", ec="gold", fa=0.5, ew=2.5)
end

for x in wit.mlist_inside[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs_)::Int
    # plot_point!(ax_[iloc], x[1:2], mc="blue")
end

for x in wit.mlist_image[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs_)::Int
    plot_point!(ax_[iloc], x[1:2], mc="purple")
end

for x in wit.mlist_outside[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs_)::Int
    plot_point!(ax_[iloc], x[1:2], mc="red")
end

for x in wit.mlist_unknown[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs_)::Int
    plot_point!(ax_[iloc], x[1:2], mc="orange")
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