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

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

N = 6
M = 1

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
τ0 = 1.5 - dt

γ = 100
TransOFF = [γ, 0, 0, 0]
TransON  = [0, γ, 0, 0]
StabOFF  = [0, 0, γ, 0]
StabON   = [0, 0, 0, γ]

ZN() = zeros(N, N)
α = exp(-cf*dt)
βlo = Tlo*(1 - α)
βup = Tup*(1 - α)
A_Trans = ZN()
A_Trans[1:2, 1:2] = [α 0; 0 1]
A_Stab = ZN()
A_Stab[1:2, 1:2]  = [α 0; 0 0]
b_TransOFF = [βlo, dt, TransOFF...]
b_TransON  = [βup, dt, TransON...]
b_StabOFF  = [βlo, τ0 + dt, StabOFF...]
b_StabON   = [βup, τ0 + dt, StabON...]

dom_above_Llo = [AffForm([-1, 0, 0, 0, 0, 0], Llo)]
dom_below_Llo = [AffForm([1, 0, 0, 0, 0, 0], -Llo)]
dom_above_Lup = [AffForm([-1, 0, 0, 0, 0, 0], Lup)]
dom_below_Lup = [AffForm([1, 0, 0, 0, 0, 0], -Lup)]
dom_Trans = [AffForm([0, 1, 0, 0, 0, 0], -τ0)]
dom_Stab  = [AffForm([0, -1, 0, 0, 0, 0], τ0)]
dom_TransOFF = [
    AffForm([0, 0, -1, 0, 0, 0], γ), AffForm([0, 0, 1, 0, 0, 0], -γ),
    AffForm([0, 0, 0, -1, 0, 0], 0), AffForm([0, 0, 0, 1, 0, 0], -0),
    AffForm([0, 0, 0, 0, -1, 0], 0), AffForm([0, 0, 0, 0, 1, 0], -0),
    AffForm([0, 0, 0, 0, 0, -1], 0), AffForm([0, 0, 0, 0, 0, 1], -0)
]
dom_TransON = [
    AffForm([0, 0, -1, 0, 0, 0], 0), AffForm([0, 0, 1, 0, 0, 0], -0),
    AffForm([0, 0, 0, -1, 0, 0], γ), AffForm([0, 0, 0, 1, 0, 0], -γ),
    AffForm([0, 0, 0, 0, -1, 0], 0), AffForm([0, 0, 0, 0, 1, 0], -0),
    AffForm([0, 0, 0, 0, 0, -1], 0), AffForm([0, 0, 0, 0, 0, 1], -0)
]
dom_StabOFF = [
    AffForm([0, 0, -1, 0, 0, 0], 0), AffForm([0, 0, 1, 0, 0, 0], -0),
    AffForm([0, 0, 0, -1, 0, 0], 0), AffForm([0, 0, 0, 1, 0, 0], -0),
    AffForm([0, 0, 0, 0, -1, 0], γ), AffForm([0, 0, 0, 0, 1, 0], -γ),
    AffForm([0, 0, 0, 0, 0, -1], 0), AffForm([0, 0, 0, 0, 0, 1], -0)
]
dom_StabON = [
    AffForm([0, 0, -1, 0, 0, 0], 0), AffForm([0, 0, 1, 0, 0, 0], -0),
    AffForm([0, 0, 0, -1, 0, 0], 0), AffForm([0, 0, 0, 1, 0, 0], -0),
    AffForm([0, 0, 0, 0, -1, 0], 0), AffForm([0, 0, 0, 0, 1, 0], -0),
    AffForm([0, 0, 0, 0, 0, -1], γ), AffForm([0, 0, 0, 0, 0, 1], -γ)
]

# To TransOFF
pf_dom = PolyFunc([dom_TransOFF..., dom_above_Llo..., dom_Trans...])
piece_TOFF2TOFF = Piece(pf_dom, 1, A_Trans, b_TransOFF, 1)
pf_dom = PolyFunc([dom_TransON..., dom_above_Lup..., dom_Trans...])
piece_TON2TOFF = Piece(pf_dom, 1, A_Trans, b_TransOFF, 1)
# To TransON
pf_dom = PolyFunc([dom_TransOFF..., dom_below_Llo..., dom_Trans...])
piece_TOFF2TON = Piece(pf_dom, 1, A_Trans, b_TransON, 1)
pf_dom = PolyFunc([dom_TransON..., dom_below_Lup..., dom_Trans...])
piece_TON2TON = Piece(pf_dom, 1, A_Trans, b_TransON, 1)
# To StabOFF
pf_dom = PolyFunc([dom_TransOFF..., dom_above_Llo..., dom_Stab...])
piece_TOFF2SOFF = Piece(pf_dom, 1, A_Stab, b_StabOFF, 1)
pf_dom = PolyFunc([dom_TransON..., dom_above_Lup..., dom_Stab...])
piece_TON2SOFF = Piece(pf_dom, 1, A_Stab, b_StabOFF, 1)
pf_dom = PolyFunc([dom_StabOFF..., dom_above_Llo...])
piece_SOFF2SOFF = Piece(pf_dom, 1, A_Stab, b_StabOFF, 1)
pf_dom = PolyFunc([dom_StabON..., dom_above_Lup...])
piece_SON2SOFF = Piece(pf_dom, 1, A_Stab, b_StabOFF, 1)
# To StabON
pf_dom = PolyFunc([dom_TransOFF..., dom_below_Llo..., dom_Stab...])
piece_TOFF2SON = Piece(pf_dom, 1, A_Stab, b_StabON, 1)
pf_dom = PolyFunc([dom_TransON..., dom_below_Lup..., dom_Stab...])
piece_TON2SON = Piece(pf_dom, 1, A_Stab, b_StabON, 1)
pf_dom = PolyFunc([dom_StabOFF..., dom_below_Llo...])
piece_SOFF2SON = Piece(pf_dom, 1, A_Stab, b_StabON, 1)
pf_dom = PolyFunc([dom_StabON..., dom_below_Lup...])
piece_SON2SON = Piece(pf_dom, 1, A_Stab, b_StabON, 1)
#
sys = System([
    piece_TOFF2TOFF, piece_TON2TOFF,
    piece_TOFF2TON, piece_TON2TON,
    piece_TOFF2SOFF, piece_TON2SOFF, piece_SOFF2SOFF, piece_SON2SOFF,
    piece_TOFF2SON, piece_TON2SON, piece_SOFF2SON, piece_SON2SON
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
x = [T0, 0.0, TransOFF...]
tol_dom = 1e-8

for istep = 1:nstep
    global x
    T = x[1]
    loc = x[3:6]
    @assert any(l -> norm(loc - l) < 1e-5, (TransOFF, TransON, StabOFF, StabON))
    color = any(l -> norm(loc - l) < 1e-5, (TransOFF, StabOFF)) ? "blue" : "red"
    marker = any(l -> norm(loc - l) < 1e-5, (TransOFF, TransON)) ? "x" : "."
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
    AffForm([-1, 0, 0, 0, (Slo - Tlo)/γ, (Slo - Tlo)/γ], Tlo),
    AffForm([1, 0, 0, 0, (Tup - Sup)/γ, (Tup - Sup)/γ], -Tup),
    AffForm([0, -1, 0, 0, 0, 0], 0),
    AffForm([0, 1, 0, 0, 0, 0], -(τ0 + 2*dt))
])])
# empty!(mpf_safe.pfs[3].afs)
# empty!(mpf_safe.pfs[4].afs)

mlist_init = [[
    [Ilo, 0, TransOFF...], [Iup, 0, TransOFF...],
    [Ilo, 0, TransON...], [Iup, 0, TransON...]
]]

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
ylims = (-5*dt, τ0 + 5*dt)
lims = ([Tlo - 2*DT, -τ0 - 5*dt], [Tup + 2*DT, 2*τ0 + 5*dt])

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.tick_params(axis="both", labelsize=15)
end

afs = AffForm{Vector{Float64},Float64}[]
locs = (TransOFF, TransON, StabOFF, StabON)

for (iloc, loc) in enumerate(locs)
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
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs)::Int
    # plot_point!(ax_[iloc], x[1:2], mc="blue")
end

for x in wit.mlist_image[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs)::Int
    plot_point!(ax_[iloc], x[1:2], mc="purple")
end

for x in wit.mlist_outside[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs)::Int
    plot_point!(ax_[iloc], x[1:2], mc="red")
end

for x in wit.mlist_unknown[1]
    loc = x[3:6]
    iloc = findfirst(l -> norm(loc - l) < 1e-5, locs)::Int
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