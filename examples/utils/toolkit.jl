using LinearAlgebra
using JuMP
using Gurobi
using PyPlot
using PyCall
const spatial = pyimport_conda("scipy.spatial", "scipy")
const optimize = pyimport_conda("scipy.optimize", "scipy")

include("../../src/CEGISPolyhedralBarrier.jl")
CPB = CEGISPolyhedralBarrier
AffForm = CPB.AffForm
GenForm = CPB.GenForm
Piece = CPB.Piece
State = CPB.State
BarrierProblem = CPB.BarrierProblem
isless_margin = CPB.isless_margin

const GUROBI_ENV = Gurobi.Env()
solver() = Model(optimizer_with_attributes(
    () -> Gurobi.Optimizer(GUROBI_ENV), "OutputFlag"=>false
))

_isin(afs, x, rtol) = all(af -> isless_margin(af, x, 0, rtol), afs)

function next_state(pieces, state, tol_dom)
    for piece in pieces
        if piece.loc_src == state.loc && _isin(piece.afs_dom, state.x, tol_dom)
            return State(piece.loc_dst, piece.A*state.x + piece.b)
        end
    end
    return State(-1, fill(NaN, length(state.x)))
end

#-------------------------------------------------------------------------------
function compute_vertices_hrep(A, b)
    @assert (size(A, 1),) == size(b)
    nvar = size(A, 2)
    M = hcat(A, -b)
    A_ub = hcat(A, map(r -> norm(r), eachrow(A)))
    c_obj = zeros(nvar + 1)
    c_obj[nvar + 1] = -1
    bounds = ((nothing, nothing), (nothing, nothing), (nothing, 1))
    res = optimize.linprog(c_obj, A_ub=A_ub, b_ub=b, bounds=bounds)
    @assert res["success"] && res["status"] == 0
    res["fun"] > 0 && return Vector{Float64}[]
    x = res["x"][1:nvar]
    hs = spatial.HalfspaceIntersection(M, x)
    points = collect.(eachrow(hs.intersections))
    ch = spatial.ConvexHull(points)
    return [ch.points[i + 1, :] for i in ch.vertices]
end

function _plot_hrep2D!(ax, A, b, fc, fa, ec, ew)
    verts = compute_vertices_hrep(A, b)
    isempty(verts) && return
    polylist = matplotlib.collections.PolyCollection((verts,))
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
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
                       fc="green", fa=0.5, ec="green", ew=1.0)
    @assert all(af -> size(af.a) == (2,), afs)
    nrow = length(afs) + 4
    A = zeros(nrow, 2)
    b = zeros(nrow)
    for (i, af) in enumerate(afs)
        A[i, :] = af.a
        b[i] = -af.β
    end
    for i = 1:2
        A[length(afs) + 1 + 2*(i - 1), i] = -1
        b[length(afs) + 1 + 2*(i - 1)] = -lims[1][i]
        A[length(afs) + 2 + 2*(i - 1), i] = +1
        b[length(afs) + 2 + 2*(i - 1)] = +lims[2][i]
    end
    _plot_hrep2D!(ax, A, b, fc, fa, ec, ew)
end

function plot_level2D!(ax, gfs::Vector{GenForm}, loc, lims;
                       fc="green", fa=0.5, ec="green", ew=1.0)
    afs = extract_afs(gfs, loc)
    plot_level2D!(ax, afs, lims, fc=fc, fa=fa, ec=ec, ew=ew)
end

function plot_point!(ax, point; mc="blue", ms=15)
    ax.plot(point..., marker=".", ms=ms, c=mc)
end