module ExampleTrainsDistance

include("../toolkit.jl")
include("./problem_distance.jl")

using Random
Random.seed!(0)
colors = collect(keys(matplotlib.colors.TABLEAU_COLORS))

Nt = 4
α = 0.1
β = 0.5
lim_up = +0.02
lim_lo = -0.02
vinit = 0.1
xsafe_lo = -1.0
xsafe_up = +1.0
vsafe = 0.3

prob = build_problem(Nt, α, β, lim_lo, lim_up, vinit,
                     xsafe_lo, xsafe_up, vsafe)

# simulation
fig = figure(0, figsize=(10, 5))
axx = fig.add_subplot()
axv = axx.twinx()

nstep = 15
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

for state_init in prob.states_init
    state = state_init
    @assert abs(sum(state.x)) < 1e-8
    x_traj = [state.x]
    for t = 1:nstep
        @assert state.loc > 0
        state = next_state(prob.pieces, state, 1e-6)
        push!(x_traj, state.x)
    end
    for i = 1:Nt
        c = colors[mod(i - 1, length(colors)) + 1]
        times = 0:nstep
        axx.plot(times, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
        axv.plot(times, getindex.(x_traj, Nt + i), ls="--", marker=".", c=c)
    end
end

# Solve !!!
iter_max = Inf
# iter_max = 5

for Nt = 4:5
    println("Nt: ", Nt)
    prob = build_problem(Nt, α, β, lim_lo, lim_up, vinit,
                         xsafe_lo, xsafe_up, vsafe)
    println("# pieces: ", length(prob.pieces))
    println("# initial states: ", length(prob.states_init))
    status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max, solver,
                                                   print_period=100)
    display(status)
    @assert Int(status) ∈ (1, 3)
    ## Algo illustration
    fig = figure(1 + Nt, figsize=(10, 5))
    ax = fig.add_subplot()
    ax2 = ax.twinx()
    ax2.set_yscale("log")

    ax.plot(rec.ninside)
    ax.plot(rec.nimage)
    ax.plot(rec.nunknown)
    ax2.plot(rec.times, ls="dashed", lw=2, c="k")
end

end # module