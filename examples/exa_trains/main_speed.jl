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

include("./problem_speed.jl")

function next_state(sys, x)
    for piece in sys.pieces
        !CPB._prox(piece.pf_dom, x, 1e-6) && continue
        return piece.A*x + piece.b
    end
    return nothing
end

Nt = 5
α = 0.5
Tstab = 12
lim_up = +0.02
lim_lo = -0.02
xinit = 0.1
xsafes = [0.3, 0.05]

N, M, sys, mlist_init, mpf_inv, mpf_safe =
    build_problem(Nt, α, Tstab, lim_lo, lim_up, xinit, xsafes)

# simulation
xsample = 0.1
nsample = 0

fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 15
ax.set_xlim((0, nstep))
ax.set_xlabel("time")

ax.set_ylim(-0.31, 0.31)
ax.plot((0, nstep), (0, 0), c="k")
ax.set_ylabel("x")

for xsafe in xsafes
    ax.plot((0, nstep), (+xsafe, +xsafe), c="r")
    ax.plot((0, nstep), (-xsafe, -xsafe), c="r")
end

for (s, x) in enumerate(mlist_init[1])
    c = colors[mod(s - 1, length(colors)) + 1]
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        push!(x_traj, x)
    end
end

for s = 1:nsample
    c = colors[mod(s - 1, length(colors)) + 1]
    xt = rand(Nt)
    xt = xt .- sum(xt)/Nt
    xt = xt/norm(xt, Inf)*xsample
    x = [xt; 0.0]
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        push!(x_traj, x)
    end
end

f = open(string(@__DIR__, "/traj_speed_", Nt, ".txt"), "w")
println(f, [string("y", i, ",") for i = 1:Nt]..., "t")
x = mlist_init[1][1]
@assert abs(sum(x)) < 1e-8
x_traj = [x]
println(f, [string(x[i], ",") for i = 1:Nt]..., "0")
for t = 1:nstep
    global x
    x = next_state(sys, x)
    isnothing(x) && error("x is nothing")
    noise = rand(N)
    x = x + (noise .- sum(noise)/N)*0.0
    push!(x_traj, x)
    println(f, [string(x[i], ",") for i = 1:Nt]..., t)
end
for i = 1:Nt
    c = colors[mod(i - 1, length(colors)) + 1]
    times = 0:nstep
    ax.plot(times, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
end
close(f)

# Solve !!!
ϵ = xsafes[2]*α/5
δ = 1e-8
iter_max = Inf

for Nt = 1:5
    N, M, sys, mlist_init, mpf_inv, mpf_safe =
        build_problem(Nt, α, Tstab, lim_lo, lim_up, xinit, xsafes)
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