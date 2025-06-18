module ExampleSaturation

# Closed loop analysis of control command software
# Equation 4
# https://hal.science/hal-01521987/file/Roux2015.pdf

include("../toolkit.jl")

# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

N = 4
A = [
    0.499 -0.05 0 0;
    0.01 1 0 0;
    0.028224 0 1 0.01;
    5.6448 0 -0.01 1;
]
b = [1, 0, -0.064, -12.8]
c = [0, 0, 1, 0]
urange = (-0.5, 0.0, 0.5)
xmaxs = (2, 1, 1, 25)
α = 1.0 # parameter [0,1] to tune: the closer to 0, the easier

gfs_safe = GenForm[]
for (i, xmax) in enumerate(xmaxs)
    push!(gfs_safe, GenForm(1, AffForm(-hot(i, N), -xmax)))
    push!(gfs_safe, GenForm(1, AffForm(+hot(i, N), -xmax)))
end

pieces = Piece[]
for u in urange
    afs_dom = [AffForm(+c, -u + 1)] # C*x - u < -1 => -1
    push!(pieces, Piece(afs_dom, 1, α*A, -α*b, 1))
    afs_dom = [AffForm(-c, +u - 1), AffForm(+c, -u - 1)] # -1 ≤ C*x - u ≤ 1
    push!(pieces, Piece(afs_dom, 1, α*A + α*b*c', -α*b*u, 1))
    afs_dom = [AffForm(-c, +u + 1)] # C*x - u > +1 => +1
    push!(pieces, Piece(afs_dom, 1, α*A, +α*b, 1))
end

# simulation
fig = figure(0, figsize=(10, 5))
ax_ = fig.subplots(
    nrows=2, ncols=2,
    gridspec_kw=Dict("wspace"=>0.2, "hspace"=>0.1)
)

nstep = 300
tol_dom = 1e-8

for (i, xmax) in enumerate(xmaxs)
    ax_[i].set_xlim((0, nstep))
    ax_[i].set_xlabel("step")
    ax_[i].set_ylim(-xmax, +xmax)
    ax_[i].plot((0, nstep), (-xmax, -xmax), c="k")
    ax_[i].plot((0, nstep), (+xmax, +xmax), c="k")
    ax_[i].set_ylabel("x[$(i)]")
end

for _ = 1:10
    xs = Vector{Float64}[]
    x = zeros(N)
    push!(xs, x)
    for istep = 1:nstep
        pieces_valid = Piece[]
        for piece in pieces
            if _isin(piece.afs_dom, x, tol_dom)
                push!(pieces_valid, piece)
            end
        end
        @assert !isempty(pieces_valid)
        piece = rand(pieces_valid)
        x = piece.A*x + piece.b
        push!(xs, x)
    end
    for i = 1:N
        ax_[i].plot(0:nstep, getindex.(xs, i))
    end
end

prob = BarrierProblem(
    N, pieces,
    GenForm[], # gfs_inv
    gfs_safe,
    [State(1, zeros(N))], # states_init
    0.01, # ϵ
    1e-8 # δ
)

iter_max = Inf
status, gen_prob, rec = @time CPB.find_barrier(prob, iter_max,
                                               solver, print_period=10)
display(status)
@assert status == CPB.BARRIER_FOUND

## Algo illustration
fig = figure(2, figsize=(10, 5))
ax = fig.add_subplot()

ax.plot(rec.ninside)
ax.plot(rec.nimage)
ax.plot(rec.nunknown)

end # module