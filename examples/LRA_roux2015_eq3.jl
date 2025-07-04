module LRA_roux2015_eq3

# Roux et al., 2015, "Closed loop analysis of control command software"
# Equation 3
# https://hal.science/hal-01521987/file/Roux2015.pdf

include("./_preamble_.jl")

hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

N = 4
A = [
    0.499 -0.05 1 0;
    0.01 1 0 0;
    0.028224 0 0.936 0.01;
    5.6448 0 -12.81 1;
]
b = [-1, 0, 0.064, 12.8]
urange = (-0.5, 0.0, 0.5)
xmaxs = (2, 1, 1, 25)
α = 0.9 # parameter [0,1] to tune: the closer to 0, the easier

display(α * A)
display(α * b)

gfs_safe = GenForm[]
for (i, xmax) in enumerate(xmaxs)
    push!(gfs_safe, GenForm(1, AffForm(-hot(i, N), -xmax)))
    push!(gfs_safe, GenForm(1, AffForm(+hot(i, N), -xmax)))
end

afs_dom = AffForm[]
pieces = Piece[]
for u in urange
    push!(pieces, Piece(afs_dom, 1, α * A, α * b * u, 1))
end

# Illustration partial
nstep = 5
ax_list = [
    plot(xlabel="time", ylabel=string("x[", i, "]"),
         xlims=[0, nstep],
         ylims=[-xmaxs[i], +xmaxs[i]])
    for i = 1:N
]
states_init = [State(1, zeros(N))]
trajectories = build_trajectories(pieces, states_init, nstep, 1e-6)
for i = 1:N
    plot_timeseries!(ax_list[i], trajectories, i)
end

display(plot(ax_list...))

# Solve !!!
prob = BarrierProblem(
    N, pieces,
    GenForm[], # gfs_inv
    gfs_safe,
    states_init,
    0.05, # ϵ
)

iter_max = Inf
status, gen_prob = @time CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

end # module