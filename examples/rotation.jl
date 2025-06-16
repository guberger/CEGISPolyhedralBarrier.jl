module ExampleRotation

include("./utils/toolkit.jl")

N = 2

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], -2.0)),
    GenForm(1, AffForm([+1.0, 0.0], -2.0)),
    GenForm(1, AffForm([0.0, -1.0], -2.0)),
    GenForm(1, AffForm([0.0, +1.0], -2.0))
]

afs_dom = [
    AffForm([0.0, +1.0], 0.0)
]
θ = π / 5
A = [cos(θ) -sin(θ); sin(θ) cos(θ)]
b = [0.0, 0.0]
piece1 = Piece(afs_dom, 1, A, b, 1)
#
afs_dom = [
    AffForm([0.0, -1.0], 0.0)
]
μ = 0.2
piece2 = Piece(afs_dom, 1, A - μ * I, b, 1)
#
pieces = [piece1, piece2]

gfs_safe = [
    GenForm(1, AffForm([-1.0, 0.0], -2)),
    GenForm(1, AffForm([+1.0, 0.0], -2)),
    GenForm(1, AffForm([0.0, -1.0], -2)),
    GenForm(1, AffForm([0.0, +1.0], -2))
]

states_init = [
    State(1, [-1.0, -1.0]), State(1, [+1.0, -1.0]),
    State(1, [-1.0, +1.0]), State(1, [+1.0, +1.0])
]

gfs_init = [
    GenForm(1, AffForm([-1.0, 0.0], -1)),
    GenForm(1, AffForm([+1.0, 0.0], -1)),
    GenForm(1, AffForm([0.0, -1.0], -1)),
    GenForm(1, AffForm([0.0, +1.0], -1))
]

# Illustration partial
ax = plot(xlabel="x", ylabel="y",
          xlims=[-2.1, +2.1], ylims=[-2.1, +2.1],
          aspect_ratio=:equal, dpi=500)
lims = [(-10, -10), (+10, +10)]
plot!(ax, [0], [0], markershape=:xcross, ms=7, c=:black)
for gf in gfs_safe
    af = gf.af
    af_out = AffForm(-af.a, -af.β)
    plot_level2D!(ax, [af_out], lims, fa=1, fc=:red, ew=0)
end
# plot_level2D!(ax, gfs_safe, 1, lims, fa=0, ec=:green)
# plot_level2D!(ax, gfs_inv, 1, lims, fa=0, ec=:yellow)
for piece in pieces
    @assert piece.loc_src == 1
    plot_level2D!(ax, piece.afs_dom, lims, fa=0, ec=:green)
end
plot_level2D!(ax, gfs_init, 1, lims, fc=:blue, fa=0.1, ec=:blue)
for state in states_init
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:circle, c=:blue, ms=5)
end
trajectories = build_trajectories(pieces, states_init, 20, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)

savefig(ax, "examples/figures/fig_rotation_partial.png")
display(ax)

# Illustration full
ax = plot(xlabel="x", ylabel="y",
          xlims=[-2.1, +2.1], ylims=[-2.1, +2.1],
          aspect_ratio=:equal, dpi=500)
lims = [(-10, -10), (+10, +10)]
plot!(ax, [0], [0], markershape=:xcross, ms=7, c=:black)
for gf in gfs_safe
    af = gf.af
    af_out = AffForm(-af.a, -af.β)
    plot_level2D!(ax, [af_out], lims, fa=1, fc=:red, ew=0)
end
for piece in pieces
    @assert piece.loc_src == 1
    plot_level2D!(ax, piece.afs_dom, lims, fa=0, ec=:green)
end

display(ax)

# Solve !!!
prob = BarrierProblem(
    N, pieces,
    gfs_inv,
    gfs_safe,
    states_init,
    0.1, # ϵ
)

iter_max = Inf
status, gen_prob = CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

plot_level2D!(ax, gen_prob.gfs, 1, lims,
              fc="gold", ec="gold", fa=0.5)

for state in gen_prob.states_inside
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:blue, ms=5)
end

for state in gen_prob.states_image
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:purple, ms=5)
end

for edge in gen_prob.edges_unknown
    state = edge.src
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:red, ms=5)
end

states_init_traj = vcat(states_init, gen_prob.states_inside)
trajectories = build_trajectories(pieces, states_init_traj, 20, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)

savefig(ax, "examples/figures/fig_rotation_full.png")
display(ax)

end # module