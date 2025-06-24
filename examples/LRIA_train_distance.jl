module TrainsDistance

include("./_preamble_.jl")
include("./_train_distance.jl")

Nt = 3
α = 0.1
β = 0.5
lim_up = +0.02
lim_lo = -0.02
vinit = 0.1
xsafe_lo = -1.0
xsafe_up = +1.0
vsafe = 0.3

prob = build_problem(Nt, α, β, lim_lo, lim_up, vinit, xsafe_lo, xsafe_up, vsafe)

# Illustration partial
nstep = 15
ax_x_list = [
    plot(xlabel="time", ylabel=string("x[", i, "]"),
         xlims=[0, nstep],
         ylims=[-1.2, +1.2])
    for i = 1:Nt
]
ax_v_list = [
    plot(xlabel="time", ylabel=string("v[", i, "]"),
         xlims=[0, nstep],
         ylims=[-0.5, +0.5])
    for i = 1:Nt
]
trajectories = build_trajectories(prob.pieces, prob.states_init, nstep, 1e-6)
for i = 1:Nt
    plot_timeseries!(ax_x_list[i], trajectories, i)
    plot_timeseries!(ax_v_list[i], trajectories, Nt + i)
end

display(plot(ax_x_list..., ax_v_list...))

# Solve !!!
iter_max = Inf
status, gen_prob = @time CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

end # module