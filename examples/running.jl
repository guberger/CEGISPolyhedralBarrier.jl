module ExampleRunning

# Counter-example guided synthesis of neural network Lyapunov functions for piecewise linear systems
# Example 1 (Section A)
# https://ieeexplore.ieee.org/document/9304201

include("./utils/toolkit.jl")

N = 2

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], -2.0)),
    GenForm(1, AffForm([1.0, 0.0], -2.0)),
    GenForm(1, AffForm([0.0, -1.0], -2.0)),
    GenForm(1, AffForm([0.0, 1.0], -2.0))
]

afs_dom = [
    AffForm([-1.0, 0.0], 0.0), AffForm([0.0, 1.0], 0.0)
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
    AffForm([1.0, 0.0], 0.0), AffForm([0.0, 1.0], 0.0)
]
A = [-0.457 0.215; 0.491 0.49]
b = [0.0, 0.0]
piece3 = Piece(afs_dom, 1, A, b, 1)
#
afs_dom = [
    AffForm([1.0, 0.0], 0.0), AffForm([0.0, -1.0], 0.0)
]
A = [-0.022 0.344; 0.458 0.271]
b = [0.0, 0.0]
piece4 = Piece(afs_dom, 1, A, b, 1)
#
pieces = [piece1, piece2, piece3, piece4]

gfs_save = [
    GenForm(1, AffForm([-1.0, -1.0], -1.8)),
    GenForm(1, AffForm([-1.0, 1.0], -1.8)),
    GenForm(1, AffForm([1.0, -1.0], -1.8)),
    GenForm(1, AffForm([1.0, 1.0], -1.8))
]

states_init = [
    State(1, [-1.0, 0.0]), State(1, [1.0, 0.0]),
    State(1, [0.0, -1.0]), State(1, [0.0, 1.0])
]

# Illustration
ax = plot(xlims=[-3.2, 3.2], ylims=[-3.2, 3.2])
lims = [(-10, -10), (10, 10)]
plot!(ax, [0], [0], markershape=:xcross, ms=7, c=:black)
plot_level2D!(ax, gfs_save, 1, lims, fc=:green, fa=0.1, ec=:green)
plot_level2D!(ax, gfs_inv, 1, lims, fc=false, fa=0, ec=:yellow)
for piece in pieces
    @assert piece.loc_src == 1
    plot_level2D!(ax, piece.afs_dom, lims, fc=:blue, fa=0.1, ec=:blue)
end
for state in states_init
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:yellow, ms=5)
end

display(ax)

# Solve !!!
prob = BarrierProblem(
    N, pieces,
    gfs_inv,
    gfs_save,
    states_init,
    0.1, # ϵ
)

iter_max = Inf
status, gen_prob = CPB.find_barrier(prob, iter_max, solver)
display(status)
@assert status == CPB.BARRIER_FOUND

plot_level2D!(ax, gen_prob.gfs, 1, lims,
              fc="gold", ec="gold", fa=0.5, ew=2.5)

for state in gen_prob.states_inside
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:blue, ms=5)
end

for state in gen_prob.states_image
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:purple, ms=5)
end

for state in gen_prob.states_outside
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:black, ms=5)
end

for edge in gen_prob.edges_unknown
    state = edge.src
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:orange, ms=5)
end

for edge in gen_prob.edges_unknown_new
    state = edge.src
    @assert state.loc == 1
    plot!(ax, [state.x[1]], [state.x[2]], marker=:dot, c=:red, ms=5)
end

@assert isempty(gen_prob.states_outside)

display(ax)

end # module