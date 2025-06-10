using LinearAlgebra
using JuMP
using Gurobi
using Plots
using Polyhedra

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
GenForm = CPB.GenForm
Piece = CPB.Piece
State = CPB.State
BarrierProblem = CPB.BarrierProblem

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

_isin(afs, x, rtol) = all(af -> CBF._eval(af, x) ≤ rtol*norm(af), afs)

function next_state(pieces, state, tol_dom)
    for piece in pieces
        if piece.loc_src == state.loc && _isin(piece.afs_dom, state.x, tol_dom)
            return State(piece.loc_dst, piece.A*state.x + piece.b)
        end
    end
    return State(-1, fill(NaN, length(state.x)))
end

#-------------------------------------------------------------------------------

function _plot_poly2D!(ax, P, fc, fa, ec, ew)
    plot!(ax, P, fc=fc, fa=fa, lc=ec, lw=ew)
end

function build_polyhedron(afs, lims, N)
    nrow = length(afs) + 2 * N
    A = zeros(nrow, N)
    b = zeros(nrow)
    for (i, af) in enumerate(afs)
        A[i, :] = af.a
        b[i] = -af.β
    end
    for i = 1:N
        A[length(afs) + 1 + 2 * (i - 1), i] = -1
        b[length(afs) + 1 + 2 * (i - 1)] = -lims[1][i]
        A[length(afs) + 2 + 2 * (i - 1), i] = +1
        b[length(afs) + 2 + 2 * (i - 1)] = +lims[2][i]
    end
    return polyhedron(hrep(A, b))
end

function extract_afs(gfs, loc)
    afs = AffForm[]
    for gf in gfs
        if gf.loc == loc
            push!(afs, gf.af)
        end
    end
    return afs
end

function plot_level2D!(ax, afs::Vector{AffForm}, lims;
                       fc=:green, fa=0.5, ec=:green, ew=1.0)
    @assert all(af -> size(af.a) == (2,), afs)
    P = build_polyhedron(afs, lims, 2)
    _plot_poly2D!(ax, P, fc, fa, ec, ew)
end

function plot_level2D!(ax, gfs::Vector{GenForm}, loc::Int, lims;
                       fc=:green, fa=0.5, ec=:green, ew=1.0)
    afs = extract_afs(gfs, loc)
    plot_level2D!(ax, afs, lims, fc=fc, fa=fa, ec=ec, ew=ew)
end