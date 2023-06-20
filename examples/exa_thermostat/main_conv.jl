module ExampleThermostat

include("../utils/toolkit.jl")

N = 6
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

γ = 100.0
hots = [[j == i ? 1.0 : 0.0 for j = 1:6] for i = 1:6]
loc_xs = [
    [[k == 2*(j - 1) + i ? γ : 0.0 for k = 1:4] for j = 1:2] for i = 1:2
]
guards_temp = [
    [AffForm(-hots[1], +Llo), AffForm(+hots[1], -Llo)],
    [AffForm(-hots[1], +Lup), AffForm(+hots[1], -Lup)]
]
guards_time = [AffForm(+hots[2], -Tstab), AffForm(-hots[2], +Tstab)]
α = exp(-cf*dt)
As_ = [[α 0; 0 1], [α 0; 0 0]]
As = [[As_[i] zeros(2, 4); zeros(4, 6)] for i = 1:2]
bs_ = [
    [[Tlo*(1 - α), dt], [Tlo*(1 - α), Tstab + dt]],
    [[Tup*(1 - α), dt], [Tup*(1 - α), Tstab + dt]]
]
bs = [[[bs_[i][j]..., loc_xs[i][j]...] for j = 1:2] for i = 1:2]

pieces = Piece[]
for (i1, j1, i2, j2) in Iterators.product(1:2, 1:2, 1:2, 1:2)
    j1 > j2 && continue
    afs_dom = [guards_temp[i1][i2], guards_time[j2]]
    loc_x = loc_xs[i1][j1]
    for i = 3:6
        push!(afs_dom, AffForm(+hots[i], -loc_x[i - 2]))
        push!(afs_dom, AffForm(-hots[i], +loc_x[i - 2]))
    end
    push!(pieces, Piece(afs_dom, 1, As[j2], bs[i2][j2], 1))
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

state = State(1, [Tlo, 0.0, loc_xs[1][1]...])
tol_dom = 1e-8

for istep = 1:nstep
    global state
    @assert state.loc == 1
    loc_x = state.x[3:6]
    @assert any(i -> any(j -> norm(loc_x - loc_xs[i][j]) < 1e-5, 1:2), 1:2)
    color = any(j -> norm(loc_x - loc_xs[1][j]) < 1e-5, 1:2) ? "blue" : "red"
    marker = any(i -> norm(loc_x - loc_xs[i][1]) < 1e-5, 1:2) ? "x" : "."
    ax.plot((istep - 1)*dt, state.x[1], marker=marker, ms=5, c=color)
    istep == nstep && break
    state = next_state(pieces, state, tol_dom)
end

inits = [[Ilo, 0.0], [Iup, 0.0]]
xlist_init = Vector{Float64}[]
for (i, k) in Iterators.product(1:2, 1:2)
    push!(xlist_init, [inits[k]..., loc_xs[i][1]...])
end

prob = BarrierProblem(
    N, pieces,
    State[], # gfs_inv
    [
        GenForm(1, AffForm([-1, 0, 0, 0, (Slo - Tlo)/γ, (Slo - Tlo)/γ], +Tlo)),
        GenForm(1, AffForm([+1, 0, 0, 0, (Tup - Sup)/γ, (Tup - Sup)/γ], -Tup)),
        GenForm(1, AffForm([0, -1, 0, 0, 0, 0], -1.0)),
        GenForm(1, AffForm([0, +1, 0, 0, 0, 0], -(Tstab + 2*dt)))
    ], # gfs_safe
    [State(1, x) for x in xlist_init], # states_init
    dt/3, # ϵ
    1e-8 # δ
)

iter_max = Inf
status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max, solver,
                                               print_period=10, βmax = 0.0)
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

loc_xs_ = [loc_xs[1]..., loc_xs[2]...]

for (loc, loc_x) in enumerate(loc_xs_)
    gfs = [
        GenForm(1, AffForm(
            gf.af.a[1:2]*1.0, gf.af.β + dot(gf.af.a[3:6], loc_x*1.0)
        ))
        for gf in prob.gfs_safe
    ]
    plot_level2D!(ax_[loc], gfs, 1, lims, fc="none", fa=0, ec="red", ew=1.5)
    gfs = [
        GenForm(1, AffForm(
            gf.af.a[1:2]*1.0, gf.af.β + dot(gf.af.a[3:6], loc_x*1.0)
        ))
        for gf in prob.gfs_inv
    ]
    plot_level2D!(ax_[loc], gfs, 1, lims, fc="none", ec="blue")
    gfs = [
        GenForm(1, AffForm(
            gf.af.a[1:2]*1.0, gf.af.β + dot(gf.af.a[3:6], loc_x*1.0)
        ))
        for gf in gen_prob.gfs
    ]
    plot_level2D!(ax_[loc], gfs, 1, lims, fc="gold", ec="gold", fa=0.5, ew=2.5)
end

find_loc(state) = findfirst(
    loc_x -> norm(loc_x - state.x[3:6]) < 1e-5, loc_xs_
)::Int

for state in gen_prob.states_inside
    loc = find_loc(state)
    plot_point!(ax_[loc], state.x[1:2], mc="blue")
end

for state in gen_prob.states_image
    loc = find_loc(state)
    plot_point!(ax_[loc], state.x[1:2], mc="purple")
end

for state in gen_prob.states_outside
    loc = find_loc(state)
    plot_point!(ax_[loc], state.x[1:2], mc="red")
end

for link in gen_prob.links_unknown
    loc = find_loc(link.src)
    plot_point!(ax_[loc], link.src.x[1:2], mc="orange")
end

## Algo illustration
fig = figure(2, figsize=(10, 5))
ax = fig.add_subplot()

ax.plot(rec.ninside)
ax.plot(rec.nimage)
ax.plot(rec.nunknown)
ax.plot(rec.noutside)

end # module