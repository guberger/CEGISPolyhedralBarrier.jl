module ExampleICE1D

include("../utils/toolkit.jl")

N = 1

gfs_inv = [
    GenForm(1, AffForm([-1.0], -10.0)),
    GenForm(1, AffForm([1.0], 10.0))
]

afs_dom = [
    AffForm([-1.0], 0.0),
    AffForm([1.0], -4.0)
]
A = [1.0;;]
b = [1.0]
pieces = [Piece(afs_dom, 1, A, b, 1)]

gfs_save = [
    GenForm(1, AffForm([1.0], -6.0))
]

states_init = [State(1, [0.0])]

# Solve !!!
prob = BarrierProblem(
    N, pieces,
    gfs_inv,
    gfs_save,
    states_init,
    0.04, # ϵ
    1e-8 # δ
)

iter_max = Inf
status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max,
                                               solver, print_period=10,
                                               rec_gen=true)
display(status)
@assert status == CPB.BARRIER_FOUND

display(gen_prob)

end # module