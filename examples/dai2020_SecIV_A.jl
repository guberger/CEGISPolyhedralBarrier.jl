module Dai2020_SecIV_A

# Dai et al., 2020, "Counter-example guided synthesis of neural network Lyapunov
# functions for piecewise linear systems"
# Example 1 (Section IV.A)
# https://ieeexplore.ieee.org/document/9304201

include("./toolkit.jl")

N = 2

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], -2.0)),
    GenForm(1, AffForm([+1.0, 0.0], -2.0)),
    GenForm(1, AffForm([0.0, -1.0], -2.0)),
    GenForm(1, AffForm([0.0, +1.0], -2.0))
]

afs_dom = [
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, +1.0], 0.0)
]
A = [-0.999 0.0; -0.139 0.341]
b = [0.0, 0.0]
piece1 = Piece(afs_dom, 1, A, b, 1)
#
afs_dom = [
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
]
A = [0.436 0.323; 0.388 -0.049]
b = [0.0, 0.0]
piece2 = Piece(afs_dom, 1, A, b, 1)
#
afs_dom = [
    AffForm([+1.0, 0.0], 0.0), AffForm([0.0, +1.0], 0.0)
]
A = [-0.457 0.215; 0.491 0.49]
b = [0.0, 0.0]
piece3 = Piece(afs_dom, 1, A, b, 1)
#
afs_dom = [
    AffForm([+1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
]
A = [-0.022 0.344; 0.458 0.271]
b = [0.0, 0.0]
piece4 = Piece(afs_dom, 1, A, b, 1)
#
pieces = [piece1, piece2, piece3, piece4]

gfs_safe = [
    GenForm(1, AffForm([-1.0, -1.0], -1.8)),
    GenForm(1, AffForm([-1.0, +1.0], -1.8)),
    GenForm(1, AffForm([+1.0, -1.0], -1.8)),
    GenForm(1, AffForm([+1.0, +1.0], -1.8))
]

states_init = [
    State(1, [-1.0, 0.0]), State(1, [+1.0, 0.0]),
    State(1, [0.0, -1.0]), State(1, [0.0, +1.0])
]

# Illustration
ax = plot(xlims=[-2.1, +2.1], ylims=[-2.1, +2.1], aspect_ratio=:equal)
lims = [(-10, -10), (+10, +10)]
plot!(ax, [0], [0], markershape=:xcross, ms=7, c=:black)
for gf in gfs_safe
    gf_out = GenForm(gf.loc, AffForm(-gf.af.a, -gf.af.β))
    plot_level2D!(ax, [gf_out], 1, lims, fa=1, fc=:red, lw=0)
end
# plot_level2D!(ax, gfs_safe, 1, lims, fa=0, lc=:green)
# plot_level2D!(ax, gfs_inv, 1, lims, fa=0, lc=:yellow)
for piece in pieces
    @assert piece.loc_src == 1
    plot_level2D!(ax, piece.afs_dom, lims, fa=0, lc=:blue)
end
for state in states_init
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:circle, c=:yellow, ms=5)
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

trajectories = build_trajectories(pieces, gen_prob.states_inside, 10, 1e-6)
plot_trajectories2D!(ax, trajectories, 1)
lens!([-1.18, -1.0], [-0.2, -0.1], inset=(1, bbox(0.65, 0.05, 0.3, 0.3)))

display(ax)

end # module