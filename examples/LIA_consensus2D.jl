module LIA_consensus2D

include("./_preamble_.jl")

N = 3 # x, y, t
vmax = +10
vmin = -10
T = vmax - vmin

mode_to_loc = Dict(["t"=>1, "s"=>2])
guards_xy = Dict([
    "x"=>[AffForm([-1, +1, 0], +1)],
    "y"=>[AffForm([+1, -1, 0], +1)],
    "="=>[AffForm([-1, +1, 0], 0), AffForm([+1, -1, 0], 0)]
])
guards_t = Dict([
    ("t", "t")=>[AffForm([0, 0, -1], 0), AffForm([0, 0, +1], -(T - 2))],
    ("t", "s")=>[AffForm([0, 0, -1], +(T - 1)), AffForm([0, 0, +1], -(T - 1))],
    ("s", "s")=>[AffForm([0, 0, -1], +T), AffForm([0, 0, +1], -T)]
])
As = Dict(["t"=>[1 0 0; 0 1 0; 0 0 1], "s"=>[1 0 0; 0 1 0; 0 0 0]])
bs = Dict([
    ("x", "t")=>[-1, 0, 1],
    ("x", "s")=>[-1, 0, T],
    ("y", "t")=>[0, -1, 1],
    ("y", "s")=>[0, -1, T],
    ("=", "t")=>[0, 0, 1],
    ("=", "s")=>[0, 0, T]
])

pieces = Piece[]
for (i, edge) in Iterators.product(("x", "y", "="),
                                   (("t", "t"), ("t", "s"), ("s", "s")))
    j1, j2 = edge
    afs_dom = [guards_xy[i]..., guards_t[edge]...]
    q1 = mode_to_loc[j1]
    q2 = mode_to_loc[j2]
    push!(pieces, Piece(afs_dom, q1, As[j2], bs[(i, j2)], q2))
end

η = 0.1
prob = BarrierProblem(
    N, pieces,
    GenForm[], # gfs_inv
    [
        [GenForm(loc, AffForm([-1, 0, 0], +vmin - η)) for loc = 1:2]...,
        [GenForm(loc, AffForm([+1, 0, 0], -vmax - η)) for loc = 1:2]...,
        [GenForm(loc, AffForm([0, -1, 0], +vmin - η)) for loc = 1:2]...,
        [GenForm(loc, AffForm([0, +1, 0], -vmax - η)) for loc = 1:2]...,
        GenForm(2, AffForm([+1, -1, 0], -η)),
        GenForm(2, AffForm([-1, +1, 0], -η))
    ], # gfs_safe
    [
        State(1, [vmin, vmin, 0]), State(1, [vmin, vmax, 0]),
        State(1, [vmax, vmin, 0]), State(1, [vmax, vmax, 0])
    ], # states_init
    0.01, # ϵ
)

iter_max = Inf
status, gen_prob = CPB.find_barrier(prob, iter_max, solver,
                                    δ=0, isint_gen=true, isint_verif=true)
display(status)
@assert status == CPB.BARRIER_FOUND
display(gen_prob.gfs)

end # module