module ExampleThermostat

include("../utils/toolkit.jl")

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

mode_to_locs = [[1, 3], [2, 4]] # cool [trans, stab], heat [trans, stab]
guards_temp = [
    [AffForm([-1.0, 0.0], Llo), AffForm([1.0, 0.0], -Llo)],
    [AffForm([-1.0, 0.0], Lup), AffForm([1.0, 0.0], -Lup)]
]
guards_time = [AffForm([0.0, 1.0], -Tstab), AffForm([0.0, -1.0], Tstab)]
α = exp(-cf*dt)
As = [[α 0; 0 1], [α 0; 0 0]]
bs = [
    [[Tlo*(1 - α), dt], [Tlo*(1 - α), Tstab + dt]],
    [[Tup*(1 - α), dt], [Tup*(1 - α), Tstab + dt]]
]

pieces = Piece[]
for (i1, j1, i2, j2) in Iterators.product(1:2, 1:2, 1:2, 1:2)
    j1 > j2 && continue
    afs_dom = [guards_temp[i1][i2], guards_time[j2]]
    push!(pieces, Piece(afs_dom,
                        mode_to_locs[i1][j1],
                        As[j2], bs[i2][j2],
                        mode_to_locs[i2][j2]))
end

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
state = State(2, [T0, 0.0])
tol_dom = 1e-8

for istep = 1:nstep
    global state
    @assert state.loc ∈ (1, 2, 3, 4)
    color = state.loc ∈ (1, 3) ? "blue" : "red"
    marker = state.loc ∈ (1, 2) ? "x" : "."
    ax.plot((istep - 1)*dt, state.x[1], marker=marker, ms=5, c=color)
    istep == nstep && break
    state = next_state(pieces, state, tol_dom)
end

prob = BarrierProblem(
    N, pieces,
    GenForm[], # gfs_inv
    [
        [GenForm(loc, AffForm([0.0, -1.0], -1.0)) for loc = 1:4]...,
        [GenForm(loc, AffForm([0.0, 1.0], -(Tstab + 2*dt))) for loc = 1:4]...,
        [GenForm(loc, AffForm([-1.0, 0.0], Tlo)) for loc in (1, 2)]...,
        [GenForm(loc, AffForm([1.0, 0.0], -Tup)) for loc in (1, 2)]...,
        [GenForm(loc, AffForm([-1.0, 0.0], Slo)) for loc in (3, 4)]...,
        [GenForm(loc, AffForm([1.0, 0.0], -Sup)) for loc in (3, 4)]...
    ], # gfs_safe
    [
        [State(loc, [Ilo, 0.0]) for loc in (1, 2)]...,
        [State(loc, [Iup, 0.0]) for loc in (1, 2)]...
    ], # states_init
    dt/3, # ϵ
    1e-8 # δ
)

iter_max = Inf
status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max,
                                               solver, print_period=10)
@assert status == CPB.BARRIER_FOUND

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

for loc = 1:4
    plot_level2D!(ax_[loc], prob.gfs_safe, loc, lims,
                  fc="none", fa=0, ec="red", ew=1.5)
end

for loc = 1:4
    plot_level2D!(ax_[loc], gen_prob.gfs, loc, lims,
                  fc="gold", ec="gold", fa=0.5, ew=2.5)
end

for state in gen_prob.states_inside
    plot_point!(ax_[state.loc], state.x, mc="blue")
end

for state in gen_prob.states_image
    plot_point!(ax_[state.loc], state.x, mc="purple")
end

for link in gen_prob.links_unknown
    state = link.src
    plot_point!(ax_[state.loc], state.x, mc="orange")
end

@assert isempty(gen_prob.states_outside)

## Algo illustration
fig = figure(2, figsize=(10, 5))
ax = fig.add_subplot()

ax.plot(rec.ninside)
ax.plot(rec.nimage)
ax.plot(rec.nunknown)

end # module