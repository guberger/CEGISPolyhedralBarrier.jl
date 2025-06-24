module TrainsSpeed

include("./_preamble_.jl")
include("./_train_speed.jl")

Nt = 3
α = 0.5
Tstab = 12
lim_up = +0.02
lim_lo = -0.02
vinit = 0.1
vsafes = [0.3, 0.05]

prob = build_problem(Nt, α, Tstab, lim_lo, lim_up, vinit, vsafes)

# Illustration partial
nstep = 15
ax_list = [
    plot(xlabel="time", ylabel=string("x[", i, "]"),
         xlims=[0, nstep],
         ylims=[-0.5, +0.5])
    for i = 1:Nt
]
trajectories = build_trajectories(prob.pieces, prob.states_init, nstep, 1e-6)
for i = 1:Nt
    plot_timeseries!(ax_list[i], trajectories, i)
end

display(plot(ax_list...))

# Solve !!!
iter_max = Inf
status, gen_prob = @time CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

end # module