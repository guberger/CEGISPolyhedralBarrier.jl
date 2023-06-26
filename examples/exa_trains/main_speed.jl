module ExampleTrainsSpeed

include("../utils/toolkit.jl")
include("./problem_speed.jl")

using Random
Random.seed!(0)
colors = collect(keys(matplotlib.colors.TABLEAU_COLORS))

Nt = 5
α = 0.5
Tstab = 12
lim_up = +0.02
lim_lo = -0.02
vinit = 0.1
vsafes = [0.3, 0.05]

prob = build_problem(Nt, α, Tstab, lim_lo, lim_up, vinit, vsafes)

# simulation
fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 15
ax.set_xlim((0, nstep))
ax.set_xlabel("time")

ax.plot((0, nstep), (0, 0), c="k")
ax.set_ylabel("v")

for vsafe in vsafes
    ax.plot((0, nstep), (+vsafe, +vsafe), c="r")
    ax.plot((0, nstep), (-vsafe, -vsafe), c="r")
end

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
        ax.plot(times, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
    end
end

# Solve !!!
iter_max = Inf
# iter_max = 5

for Nt = 1:1
    prob = build_problem(Nt, α, Tstab, lim_lo, lim_up, vinit, vsafes)
    println("# pieces: ", length(prob.pieces))
    println("# initial states: ", length(prob.states_init))
    status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max, solver,
                                                   print_period=10)
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
    ax.plot(rec.noutside)
    ax2.plot(rec.times, ls="dashed", lw=2, c="k")
end

end # module