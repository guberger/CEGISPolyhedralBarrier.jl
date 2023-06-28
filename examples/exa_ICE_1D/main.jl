module ExampleICE1D

include("../utils/toolkit.jl")

N = 2

gfs_inv = [
    GenForm(1, AffForm([-1.0, 0.0], 0.0)),
    GenForm(1, AffForm([0.0, -1.0], 0.0)),
    GenForm(1, AffForm([0.0, 1.0], 0.0))
]

afs_dom = [AffForm([1.0, 0.0], -1.0)]
A = [1.0 0.0; 0.0 0.0]
b = [0.01, 0.0]
pieces = [Piece(afs_dom, 1, A, b, 1)]

gfs_save = [
    GenForm(1, AffForm([1.0, 0.0], -1.1))
]

states_init = [State(1, [0.0, 0.0])]

# Illustration
fig = figure(0, figsize=(10, 5))
ax_ = fig.subplots(
    nrows=4, ncols=4,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1)
)

xlims = (-0.1, 1.2)
ylims = (-1.1, 1.1)
lims = [(-1000, -1000), (1000, 1000)]

for ax in ax_
    ax.set_xlim(xlims...)
    ax.set_ylim(ylims...)
    ax.plot(0, 0, marker="x", ms=7, c="black", mew=1.5)
    plot_level2D!(ax, gfs_save, 1, lims, fc="green", fa=0.1, ec="green")
    plot_level2D!(ax, gfs_inv, 1, lims, fc="none", ec="yellow")
    for piece in pieces
        @assert piece.loc_src == 1
        plot_level2D!(ax, piece.afs_dom, lims, fc="blue", fa=0.1, ec="blue")
    end
end

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

fig = figure(2, figsize=(10, 5))
ax = fig.add_subplot()

ax.plot(rec.ninside)
ax.plot(rec.nimage)
ax.plot(rec.nunknown)

for (iter, gen_prob) in enumerate(rec.gen_probs)
    plot_level2D!(ax_[iter], gen_prob.gfs, 1, lims,
                  fc="gold", ec="gold", fa=0.5, ew=2.5)

    for state in gen_prob.states_inside
        @assert state.loc == 1
        plot_point!(ax_[iter], state.x, mc="blue")
    end
    
    for state in gen_prob.states_image
        @assert state.loc == 1
        plot_point!(ax_[iter], state.x, mc="purple")
    end
    
    for link in gen_prob.links_unknown
        state = link.src
        @assert state.loc == 1
        plot_point!(ax_[iter], state.x, mc="orange")
    end

    for link in gen_prob.links_unknown_new
        state = link.src
        @assert state.loc == 1
        plot_point!(ax_[iter], state.x, mc="red")
    end

    @assert isempty(gen_prob.states_outside)
end

end # module