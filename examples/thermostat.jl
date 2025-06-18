module Thermostat

include("./_preamble_.jl")

N = 2
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

mode_to_loc = Dict(["dt"=>1, "ut"=>2, "ds"=>3, "us"=>4])
guards_temp = Dict([
    "dd"=>AffForm([-1.0, 0.0], +Llo),
    "du"=>AffForm([+1.0, 0.0], -Llo),
    "ud"=>AffForm([-1.0, 0.0], +Lup),
    "uu"=>AffForm([+1.0, 0.0], -Lup)
])
guards_time = Dict([
    "tt"=>AffForm([0.0, +1.0], -Tstab),
    "ts"=>AffForm([0.0, -1.0], +Tstab),
    "ss"=>AffForm([0.0, -1.0], +Tstab)
])
α = exp(-cf*dt)
As = Dict("t"=>[α 0; 0 1], "s"=>[α 0; 0 0])
bs = Dict([
    "dt"=>[Tlo*(1 - α), dt],
    "ds"=>[Tlo*(1 - α), Tstab + dt],
    "ut"=>[Tup*(1 - α), dt],
    "us"=>[Tup*(1 - α), Tstab + dt]
])

pieces = Piece[]
for (i1, i2, t_edge) in Iterators.product(
        ("d", "u"), ("d", "u"), (("t", "t"), ("t", "s"), ("s", "s")))
    j1, j2 = t_edge
    afs_dom = [guards_temp[string(i1, i2)], guards_time[string(j1, j2)]]
    q1 = mode_to_loc[string(i1, j1)]
    q2 = mode_to_loc[string(i2, j2)]
    push!(pieces, Piece(afs_dom, q1, As[j2], bs[string(i2, j2)], q2))
end

# Illustration partial
nstep = 300
DT = Tup - Tlo
ax = plot(xlabel="time", ylabel="temp",
          xlims=[0, Tstab + 3 * dt],
          ylims=[Tlo - 0.05*DT, Tup + 0.05*DT],
          dpi=500)
for y in (Tlo, Slo, Llo, Lup, Sup, Tup)
    plot!(ax, [0, nstep], [y, y], c=:black, legend=false)
end
states_init_traj = [State(2, [Tlo, 0.0]), State(2, [Tup, 0.0])]
trajectories = build_trajectories(pieces, states_init_traj, nstep, 1e-6)
plot_trajectories2D!(ax, trajectories, (2, 1), (1, 2), mc=:black)
plot_trajectories2D!(ax, trajectories, (2, 1), (3, 4), mc=:red, msc=:red)

display(ax)

# Solve !!!
Tsup = Tstab + 2 * dt
prob = BarrierProblem(
    N, pieces,
    GenForm[], # gfs_inv
    [
        [GenForm(loc, AffForm([0.0, -1.0], -1.0)) for loc = 1:4]...,
        [GenForm(loc, AffForm([0.0, +1.0], -Tsup)) for loc = 1:4]...,
        [GenForm(loc, AffForm([-1.0, 0.0], +Tlo)) for loc in (1, 2)]...,
        [GenForm(loc, AffForm([+1.0, 0.0], -Tup)) for loc in (1, 2)]...,
        [GenForm(loc, AffForm([-1.0, 0.0], +Slo)) for loc in (3, 4)]...,
        [GenForm(loc, AffForm([+1.0, 0.0], -Sup)) for loc in (3, 4)]...
    ], # gfs_safe
    [
        [State(loc, [Ilo, 0.0]) for loc in (1, 2)]...,
        [State(loc, [Iup, 0.0]) for loc in (1, 2)]...
    ], # states_init
    dt/3, # ϵ
)

iter_max = Inf
status, gen_prob = CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

# Illustration
ax_list = [
    plot(xlabel="temp", ylabel="time",
         xlims=[Tlo - 0.05*DT, Tup + 0.05*DT],
         ylims=[-5 * dt, Tstab + 5 * dt],
         title=string(loc))
    for loc = 1:4
]
lims = ([Tlo - 2*DT, -Tstab - 5*dt], [Tup + 2*DT, 2*Tstab + 5*dt])
for loc = 1:4
    plot_level2D!(ax_list[loc], prob.gfs_safe, loc, lims, fa=0, lc=:red, lw=1.5)
end

for loc = 1:4
    plot_level2D!(ax_list[loc], gen_prob.gfs, loc, lims,
                  fc=:gold, lc=:gold, fa=0.5)
end

for state in prob.states_init
    plot!(ax_list[state.loc], [state.x[1]], [state.x[2]],
          marker=:circle, c=:blue, ms=5)
end

for state in gen_prob.states_inside
    plot!(ax_list[state.loc], [state.x[1]], [state.x[2]],
          marker=:dot, c=:blue, ms=5)
end

for state in gen_prob.states_image
    plot!(ax_list[state.loc], [state.x[1]], [state.x[2]],
          marker=:dot, c=:purple, ms=5)
end

for edge in gen_prob.edges_unknown
    state = edge.src
    plot!(ax_list[state.loc], [state.x[1]], [state.x[2]],
          marker=:dot, c=:red, ms=5)
end

trajectories = build_trajectories(pieces, gen_prob.states_inside, 20, 1e-6)
for loc = 1:4
    plot_trajectories2D!(ax_list[loc], trajectories, loc)
end

display(plot(ax_list...))

end # module