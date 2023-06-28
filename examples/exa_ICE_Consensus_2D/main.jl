module ExampleICEConsensus2D

include("../utils/toolkit.jl")

N = 3 # x, y, t
vmax = 3
vmin = -3
T = 8

guards_xy = [
    [AffForm([-1, 1, 0], 1)],
    [AffForm([1, -1, 0], 1)],
    [AffForm([-1, 1, 0], 0), AffForm([1, -1, 0], 0)]
]
guards_t = [
    [[AffForm([0, 0, 1], -(T - 2))], [AffForm([0, 0, -1], T - 1)]],
    [AffForm[], AffForm[]]
]
As = [[1 0 0; 0 1 0; 0 0 1], [1 0 0; 0 1 0; 0 0 0]]
bs = [
    [[-1, 0, 1], [-1, 0, T]],
    [[0, -1, 1], [0, -1, T]],
    [[0, 0, 1], [0, 0, T]]
]

pieces = Piece[]
for (i, j1, j2) in Iterators.product(1:3, 1:2, 1:2)
    j1 > j2 && continue
    afs_dom = [guards_xy[i]..., guards_t[j1][j2]...]
    push!(pieces, Piece(afs_dom, j1, As[j2], bs[i][j2], j2))
end

# simulation
fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 10
ax.set_xlim((0, nstep))
ax.set_xlabel("time")

dv = vmax - vmin
ax.set_ylim((vmin - 0.05*dv, vmax + 0.05*dv))
for y in (vmin, vmax)
    ax.plot((0, nstep), (y, y), c="k")
end
ax.set_ylabel("x, y")

state = State(1, [vmin, vmax, 0])
tol_dom = 1e-8

for istep = 1:nstep
    global state
    @assert state.loc ∈ (1, 2)
    @assert state.loc != 2 || state.x[3] ≈ T
    marker = state.loc == 1 ? "x" : "."
    ax.plot((istep - 1), state.x[1], marker=marker, ms=5, c="blue")
    ax.plot((istep - 1), state.x[2], marker=marker, ms=5, c="red")
    istep == nstep && break
    state = next_state(pieces, state, tol_dom)
end

prob = BarrierProblem(
    N, pieces,
    [
        GenForm(1, AffForm([0, 0, -1], 0)),
        GenForm(1, AffForm([0, 0, 1], -(T - 1))),
        GenForm(2, AffForm([0, 0, 1], -T)),
        GenForm(2, AffForm([0, 0, -1], T))
    ], # gfs_inv
    [
        [GenForm(loc, AffForm([-1, 0, 0], vmin - 1)) for loc = 1:2]...,
        [GenForm(loc, AffForm([1, 0, 0], -vmax - 1)) for loc = 1:2]...,
        [GenForm(loc, AffForm([0, -1, 0], vmin - 1)) for loc = 1:2]...,
        [GenForm(loc, AffForm([0, 1, 0], -vmax - 1)) for loc = 1:2]...,
        GenForm(2, AffForm([1, -1, 0], -0.5)),
        GenForm(2, AffForm([-1, 1, 0], -0.5))
    ], # gfs_safe
    [
        State(1, [vmin, vmin, 0]), State(1, [vmin, vmax, 0]),
        State(1, [vmax, vmin, 0]), State(1, [vmax, vmax, 0])
    ], # states_init
    0.05, # ϵ
    0.0 # δ
)

iter_max = Inf
status, gen_prob, rec = CPB.find_barrier(prob, iter_max,
                                         solver, print_period=10)
@assert status == CPB.BARRIER_FOUND

# Illustration
fig = figure(1, figsize=(15, 8))
ax_ = fig.subplots(
    nrows=1, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1),
    subplot_kw=Dict("projection"=>"3d")
)

xylims = (vmin - 2, vmax + 2)
zlims = (-1, T + 1)
lims = ([-100 for i = 1:N], [100 for i = 1:N])

for ax in ax_
    ax.set_xlim(xylims...)
    ax.set_ylim(xylims...)
    ax.set_zlim(zlims...)
    ax.tick_params(axis="both", labelsize=15)
end

gfs = vcat(gen_prob.gfs, prob.gfs_inv)

for loc = 1:2
    plot_level3D!(ax_[loc], gfs, loc, lims,
                  fc="red", fa=0.1, ec="red", ew=0.5)
end

# for loc = 1:4
#     plot_level2D!(ax_[loc], gen_prob.gfs, loc, lims,
#                   fc="gold", ec="gold", fa=0.5, ew=2.5)
# end

# for state in gen_prob.states_inside
#     plot_point!(ax_[state.loc], state.x, mc="blue")
# end

# for state in gen_prob.states_image
#     plot_point!(ax_[state.loc], state.x, mc="purple")
# end

# for link in gen_prob.links_unknown
#     state = link.src
#     plot_point!(ax_[state.loc], state.x, mc="orange")
# end

# @assert isempty(gen_prob.states_outside)

## Algo illustration
fig = figure(2, figsize=(10, 5))
ax = fig.add_subplot()

ax.plot(rec.ninside)
ax.plot(rec.nimage)
ax.plot(rec.nunknown)

end # module