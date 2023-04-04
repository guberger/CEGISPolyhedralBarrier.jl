module ExampleTrains

using Random
using LinearAlgebra
using Combinatorics
using JuMP
using Gurobi
using PyPlot

Random.seed!(0)
colors = collect(keys(matplotlib.colors.TABLEAU_COLORS))

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
Piece = CPB.Piece
System = CPB.System
Witness = CPB.Witness
const _AT = AffForm{Vector{Float64},Float64}
const _PT = PolyFunc{_AT}

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

# diff(i, j, N)'*x[1:N] = x[j] - x[i]
diff(i, j, N) = [k == j ? 1.0 : k == i ? -1.0 : 0.0 for k = 1:N]
# hot(i, N)'*x[1:N] = x[i]
hot(i, N) = [j == i ? 1.0 : 0.0 for j = 1:N]

N = 6
M = 1

α = 0.75
lim = 0.1
pieces = Piece{_PT,Matrix{Float64},Vector{Float64}}[]

for cases in Iterators.product([(-1, 0, 1) for i = 1:N]...)
    afs = [AffForm(ones(N), 0.0), AffForm(-ones(N), 0.0)]
    # afs = _AT[]
    A = Matrix{Float64}(I, N, N)
    b = zeros(N)
    for (i, case) in enumerate(cases)
        iprev = mod(i - 1, N) + 1
        inext = mod(i - 0, N) + 1
        if case == -1 # xnext - xprev ≤ -lim --> -lim*α
            push!(afs, AffForm(+diff(iprev, inext, N), +lim))
            b[iprev] -= lim*α
            b[inext] += lim*α
        elseif case == 1 # xnext - xprev ≥ lim --> +lim*α
            push!(afs, AffForm(-diff(iprev, inext, N), +lim))
            b[iprev] += lim*α
            b[inext] -= lim*α
        else # -lim ≤ xnext - xprev ≤ lim --> +(xnext - xprev)*α
            push!(afs, AffForm(-diff(iprev, inext, N), -lim))
            push!(afs, AffForm(+diff(iprev, inext, N), -lim))
            A[iprev, iprev] -= α
            A[iprev, inext] += α
            A[inext, iprev] += α
            A[inext, inext] -= α
        end
    end
    if all(iszero, cases)
        # display(afs)
        display(A)
        display(b)
        display(eigvals(A))
    end
    pf_dom = PolyFunc(afs)
    push!(pieces, Piece(pf_dom, 1, A, b, 1))        
end

display(length(pieces))
sys = System(pieces)

function next_state(sys, x)
    for piece in sys.pieces
        !CPB._prox(piece.pf_dom, x, 1e-6) && continue
        return piece.A*x + piece.b
    end
    return nothing
end

# simulation
fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 15
ax.set_xlim((0, nstep))
ax.set_xlabel("time")

ax.set_ylim(-1, 1)
ax.plot((0, nstep), (0, 0), c="k")
ax.set_ylabel("x")

xinit = 0.1
Nd = N ÷ 2
x_list = Vector{Float64}[]
for sets in combinations(1:N, N - Nd)
    x = fill(float(xinit), N)
    for i in sets
        x[i] = -xinit
    end
    if N == 2*Nd
        push!(x_list, x)
    else
        for i in sets
            y = copy(x)
            y[i] = 0
            push!(x_list, y)
        end
    end
end
mlist_init = [x_list]
display(x_list)
display(length(x_list))

for (s, x) in enumerate(x_list)
    c = colors[mod(s - 1, length(colors)) + 1]
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        noise = rand(N)
        x = x + (noise .- sum(noise)/N)*0.0
        push!(x_traj, x)
    end
    for i = 1:N
        ax.plot(0:nstep, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
    end
end

for s = 1:100
    c = colors[mod(s - 1, length(colors)) + 1]
    x = rand(N)
    x = x .- sum(x)/N
    x = x/norm(x, Inf)*xinit
    @assert abs(sum(x)) < 1e-8
    x_traj = [x]
    for t = 1:nstep
        x = next_state(sys, x)
        isnothing(x) && error("x is nothing")
        noise = rand(N)
        x = x + (noise .- sum(noise)/N)*0.0
        push!(x_traj, x)
    end
    for i = 1:N
        ax.plot(0:nstep, getindex.(x_traj, i), marker=".", ms=10, c=c, lw=2)
    end
end

# ------------------------------------------------------------------------------

mpf_inv = MultiPolyFunc([PolyFunc(_AT[])])

xsafe = 0.2
afs = _AT[]
for i = 1:N
    push!(afs, AffForm(+hot(i, N), -xsafe))
    push!(afs, AffForm(-hot(i, N), -xsafe))
end
mpf_safe = MultiPolyFunc([PolyFunc(afs)])

ax.plot((0, nstep), (+xsafe, +xsafe), c="k")
ax.plot((0, nstep), (-xsafe, -xsafe), c="k")

# error("Stop here please")

ϵ = lim/5
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver, do_print=true
)

display(status)
display(wit.mlist_inside)
display(wit.mlist_image)
display(wit.mlist_unknown)
display(wit.mlist_outside)

display(next_state(sys, wit.mlist_outside[1][1]))

end # module