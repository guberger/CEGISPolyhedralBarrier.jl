module Liu2022_Fig1

# Liu et al., 2022, "Linear disjunctive invariant generation with Farkas' lemma"
# Figure 1
# https://jhc.sjtu.edu.cn/~hongfeifu/manuscripta.pdf

include("./utils/toolkit.jl")

N = 2
sm = 200

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], -sm)),
    GenForm(1, AffForm([+1.0, 0.0], -sm)),
    GenForm(1, AffForm([0.0, -1.0], -sm)),
    GenForm(1, AffForm([0.0, +1.0], -sm))
]

afs_dom = [AffForm([+1.0, 0.0], -99.0), AffForm([+1.0, 0.0], -50.0)]
A = [1.0 0.0; 0.0 1.0]
b = [1.0, 0.0]
piece1 = Piece(afs_dom, 1, A, b, 1)
afs_dom = [AffForm([+1.0, 0.0], -99.0), AffForm([-1.0, 0.0], +51.0)]
A = [1.0 0.0; 0.0 1.0]
b = [1.0, 1.0]
piece2 = Piece(afs_dom, 1, A, b, 1)
pieces = [piece1, piece2]

gfs_safe = [
    GenForm(1, AffForm([+1.0, 0.0], -101.0)),
    GenForm(1, AffForm([0.0, +1.0], -151.0))
]

states_init = [State(1, [0.0, 0.0]), State(1, [0.0, 50.0])]

# Illustration partial
ax = plot(xlabel="x", ylabel="y",
          xlims=[-0.1*sm, +1.1*sm], ylims=[-0.1*sm, +1.1*sm],
          aspect_ratio=:equal, dpi=500)
lims = [(-2*sm, -2*sm), (+2*sm, +2*sm)]
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
for state in states_init
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:circle, c=:blue, ms=5)
end
trajectories = build_trajectories(pieces, states_init, 100, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)

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
status, gen_prob = CPB.find_barrier(prob, iter_max, solver,
                                    δ=0, isint_gen=true, isint_verif=true)
display(status)
@assert status == CPB.BARRIER_FOUND
display(gen_prob.gfs)

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

savefig(ax, "examples/figures/fig_liu2022_Fig1.png")
display(ax)

end # module