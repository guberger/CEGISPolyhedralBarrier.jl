module ExampleIllustrative

using LinearAlgebra
using JuMP
using Gurobi
using PyPlot

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
PolyFunc = CPB.PolyFunc
MultiPolyFunc = CPB.MultiPolyFunc
Piece = CPB.Piece
System = CPB.System
Witness = CPB.Witness

include("../utils/plotting2D.jl")

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

N = 1
M = 2

Tlo = 0
Slo = 0.05
Llo = 0.2
Lup = 0.8
Sup = 0.95
Tup = 1
cf = 1
dt = 0.1

# cooling
pf_dom = PolyFunc([AffForm([-1.0], Llo)])
A = [1.0 - cf*dt;;]
b = [cf*dt*Tlo]
piece1 = Piece(pf_dom, 1, A, b, 1)
# start heating
pf_dom = PolyFunc([AffForm([1.0], -Llo)])
A = [1.0 - cf*dt;;]
b = [cf*dt*Tup]
piece2 = Piece(pf_dom, 1, A, b, 2)
# heating
pf_dom = PolyFunc([AffForm([1.0], -Lup)])
A = [1.0 - cf*dt;;]
b = [cf*dt*Tup]
piece3 = Piece(pf_dom, 2, A, b, 2)
# start cooling
pf_dom = PolyFunc([AffForm([-1.0], Lup)])
A = [1.0 - cf*dt;;]
b = [cf*dt*Tlo]
piece4 = Piece(pf_dom, 2, A, b, 1)
#
sys = System([piece1, piece2, piece3, piece4])

# simulation
fig = figure(0, figsize=(10, 5))
ax = fig.add_subplot()

nstep = 100
ax.set_xlim((-1, nstep + 1))

ax.set_ylim((-0.1, 1.1))
for y in (Tlo, Slo, Llo, Lup, Sup, Tup)
    ax.plot((0, nstep), (y, y), c="k")
end

T0 = (Lup + Sup)/2
x = [T0]
loc = 1
tol_dom = 1e-8

for istep = 1:nstep
    global x, loc
    @assert loc ∈ (1, 2)
    color = loc == 1 ? "red" : "blue"
    ax.plot(istep, x[1], marker=".", ms=5, c=color)
    istep == nstep && break
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !CPB._prox(piece.pf_dom, x, 0) && continue
        loc = piece.loc2
        x = piece.A*x + piece.b
        break
    end
end

mpf_inv = MultiPolyFunc([
    PolyFunc(AffForm{Vector{Float64},Float64}[]) for loc = 1:2
])

mpf_safe = MultiPolyFunc([PolyFunc([
    AffForm([-1.0], Slo),
    AffForm([1.0], -Sup),
]) for loc = 1:2])

mlist_init = [[[T0]], Vector{Float64}[]]

ϵ = 0.1
δ = 1e-8
iter_max = Inf

status, mpf, wit = CPB.learn_lyapunov!(
    sys, mpf_safe, mpf_inv, mlist_init, ϵ, δ, iter_max,
    M, N, solver, solver
)

display(status)

end # module