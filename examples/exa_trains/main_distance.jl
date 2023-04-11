module ExampleTrains

using Random
using LinearAlgebra
using Combinatorics
using JuMP
using Gurobi
using PyPlot

Random.seed!(0)
colors = collect(keys(matplotlib.colors.TABLEAU_COLORS))

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

include("./problem_distance.jl")

function next_state(sys, x)
    for piece in sys.pieces
        !CPB._prox(piece.pf_dom, x, 1e-6) && continue
        return piece.A*x + piece.b
    end
    return nothing
end

Nt = 3
α = 0.5
xmin = -0.2
lim_up = +0.02
lim_lo = -0.02
vinit = 0.1
xsafe_lo = -1.0
xsafe_up = +1.0
vsafe = 0.3

N, M, sys, mlist_init, mpf_inv, mpf_safe = build_problem(
    Nt, α, xmin, lim_lo, lim_up, vinit, xsafe_lo, xsafe_up, vsafe
)

# simulation
xsample = 0.0
vsample = 0.1
nsample = 0

fig = figure(0, figsize=(10, 5))
axx = fig.add_subplot()
axv = axx.twinx()

nstep = 30
axx.set_xlim((0, nstep))
axx.set_xlabel("time")

axx.plot((0, nstep), (0, 0), c="k")
axv.plot((0, nstep), (0, 0), c="k", ls="--")
axx.set_ylabel("x")
axv.set_ylabel("v")

axx.plot((0, nstep), (xsafe_lo, xsafe_lo), c="r")
axx.plot((0, nstep), (xsafe_up, xsafe_up), c="r")
axv.plot((0, nstep), (+vsafe, +vsafe), c="b", ls="--")
axv.plot((0, nstep), (-vsafe, -vsafe), c="b", ls="--")

for (s, x) in enumerate(mlist_init[1])
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        push!(x_traj, x)
    end
end

for s = 1:nsample
    xt = rand(Nt)
    xt = xt .- sum(xt)/Nt
    xt = xt/norm(xt, Inf)*xsample
    vt = rand(Nt)
    vt = vt .- sum(vt)/Nt
    vt = vt/norm(vt, Inf)*vsample
    x = [xt; vt]
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        push!(x_traj, x)
    end
    for i = 1:Nt
        c = colors[mod(i - 1, length(colors)) + 1]
        times = 0:nstep
        axx.plot(times, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
        axv.plot(times, getindex.(x_traj, Nt + i), marker="d", ms=10, c=c, lw=2)
    end
end

f = open(string(@__DIR__, "/traj_distance_", Nt, ".txt"), "w")
println(
    f, [string("x", i, ",") for i = 1:Nt]...,
    [string("v", i, ",") for i = 1:Nt]..., "t"
)
x = mlist_init[1][1]
display(x)
@assert abs(sum(x)) < 1e-8
x_traj = [x]
println(f, [string(x[i], ",") for i = 1:N]..., "0")
for t = 1:nstep
    global x
    x = next_state(sys, x)
    isnothing(x) && error("x is nothing")
    push!(x_traj, x)
    println(f, [string(x[i], ",") for i = 1:N]..., t)
end
close(f)
for i = 1:Nt
    c = colors[mod(i - 1, length(colors)) + 1]
    times = 0:nstep
    axx.plot(times, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
    axv.plot(times, getindex.(x_traj, Nt + i), marker="d", ms=10, c=c, lw=2)
end

# Solve !!!

ϵ = vsafe*α/5
δ = 1e-8
iter_max = Inf
# iter_max = 5

for Nt = 1:3
    N, M, sys, mlist_init, mpf_inv, mpf_safe = build_problem(
        Nt, α, xmin, lim_lo, lim_up, vinit, xsafe_lo, xsafe_up, vsafe
    )
    println("# pieces: ", length(sys.pieces))
    println("# initial states: ", sum(length, mlist_init))
    status, mpf, wit = @time CPB.learn_lyapunov!(
        sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
        M, N, solver, solver, do_print=false
    )
    display(status)
    @assert Int(status) ∈ (1, 3)
end

end # module