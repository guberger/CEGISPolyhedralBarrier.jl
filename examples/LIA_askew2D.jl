module LIA_askew2D

include("./toolkit.jl")

N = 2
sm = 100

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], -sm)),
    GenForm(1, AffForm([+1.0, 0.0], -sm)),
    GenForm(1, AffForm([0.0, -1.0], -sm)),
    GenForm(1, AffForm([0.0, +1.0], -sm))
]

afs_dom = [
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
]
A = [1.0 0.0; 0.0 1.0]
b = [1.0, 2.0]
pieces = [Piece(afs_dom, 1, A, b, 1)]

gfs_safe = [
    GenForm(1, AffForm([-2.0, 1.0], 0.0))
]

states_init = [
    State(1, [10.0, 7.0]),
    State(1, [11.0, 7.0]),
    State(1, [10.0, 8.0])
]

# Illustration partial
ax = plot(xlabel="x", ylabel="y",
          xlims=[-0.1*sm, +1.1*sm], ylims=[-0.1*sm, +1.1*sm],
          aspect_ratio=:equal, dpi=500)
lims = [(-2*sm, -2*sm), (+2*sm, +2*sm)]
plot!(ax, [0], [0], markershape=:xcross, ms=7, c=:black)
for gf in gfs_safe
    gf_out = GenForm(gf.loc, AffForm(-gf.af.a, -gf.af.β))
    plot_level2D!(ax, [gf_out], 1, lims, fa=1, fc=:red, lw=0)
end
# plot_level2D!(ax, gfs_safe, 1, lims, fa=0, lc=:green)
# plot_level2D!(ax, gfs_inv, 1, lims, fa=0, lc=:yellow)
for piece in pieces
    @assert piece.loc_src == 1
    plot_level2D!(ax, piece.afs_dom, lims, fa=0, lc=:green)
end
for state in states_init
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:circle, c=:blue, ms=5)
end
trajectories = build_trajectories(pieces, states_init, 20, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)

display(ax)

# Solve !!!
prob = BarrierProblem(
    N, pieces,
    gfs_inv,
    gfs_safe,
    states_init,
    0.0, # ϵ
)

iter_max = Inf
status, gen_prob = CPB.find_barrier(prob, iter_max, solver,
                                    δ=0, isint_gen=true, isint_verif=true)
display(status)
@assert status == CPB.BARRIER_FOUND
display(gen_prob.gfs)

plot_level2D!(ax, gen_prob.gfs, 1, lims, fc=:gold, lc=:gold, fa=0.5)

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

trajectories = build_trajectories(pieces, gen_prob.states_inside, 20, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)

display(ax)

end # module