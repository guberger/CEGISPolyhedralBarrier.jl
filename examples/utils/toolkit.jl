using LinearAlgebra
using JuMP
using Gurobi
using PyPlot
using PyCall
const spatial = pyimport_conda("scipy.spatial", "scipy")
const optimize = pyimport_conda("scipy.optimize", "scipy")
const art3d = PyObject(PyPlot.art3D)

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
    N = size(A, 2)
    A_b = hcat(A, -b)
    A_ub = hcat(A, map(r -> norm(r), eachrow(A)))
    c_obj = zeros(N + 1)
    c_obj[N + 1] = -1
    bounds = [(nothing, nothing) for i = 1:N+1]
    res = optimize.linprog(c_obj, A_ub=A_ub, b_ub=b, bounds=bounds)
    @assert res["success"] && res["status"] == 0
    res["fun"] > 0 && return Vector{Float64}[]
    x = res["x"][1:N]
    hs = spatial.HalfspaceIntersection(A_b, x)
    points = collect.(eachrow(hs.intersections))
    return spatial.ConvexHull(points)
end

function _plot_hrep2D!(ax, A, b, fc, fa, ec, ew)
    ch = compute_vertices_hrep(A, b)
    verts = [ch.points[i + 1, :] for i in ch.vertices]
    polylist = matplotlib.collections.PolyCollection((verts,))
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function _plot_hrep3D!(ax, A, b, fc, fa, ec, ew)
    ch = compute_vertices_hrep(A, b)
    verts_list = [
        [ch.points[i + 1, :] for i in simplex]
        for simplex in ch.simplices
    ]
    polylist = art3d.Poly3DCollection(verts_list)
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function build_hrep_matrices(afs, lims, N)
    nrow = length(afs) + 2*N
    A = zeros(nrow, N)
    b = zeros(nrow)
    for (i, af) in enumerate(afs)
        A[i, :] = af.a
        b[i] = -af.β
    end
    for i = 1:N
        A[length(afs) + 1 + 2*(i - 1), i] = -1
        b[length(afs) + 1 + 2*(i - 1)] = -lims[1][i]
        A[length(afs) + 2 + 2*(i - 1), i] = +1
        b[length(afs) + 2 + 2*(i - 1)] = +lims[2][i]
    end
    return A, b
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
    A, b = build_hrep_matrices(afs, lims, 2)
    _plot_hrep2D!(ax, A, b, fc, fa, ec, ew)
end

function plot_level2D!(ax, gfs::Vector{GenForm}, loc::Int, lims;
                       fc="green", fa=0.5, ec="green", ew=1.0)
    afs = extract_afs(gfs, loc)
    plot_level2D!(ax, afs, lims, fc=fc, fa=fa, ec=ec, ew=ew)
end

function plot_level3D!(ax, afs::Vector{AffForm}, lims;
                       fc="green", fa=0.5, ec="green", ew=1.0)
    @assert all(af -> size(af.a) == (3,), afs)
    A, b = build_hrep_matrices(afs, lims, 3)
    _plot_hrep3D!(ax, A, b, fc, fa, ec, ew)
end

function plot_level3D!(ax, gfs::Vector{GenForm}, loc::Int, lims;
                       fc="green", fa=0.5, ec="green", ew=1.0)
    afs = extract_afs(gfs, loc)
    plot_level3D!(ax, afs, lims, fc=fc, fa=fa, ec=ec, ew=ew)
end

function plot_point!(ax, point; mc="blue", ms=15)
    ax.plot(point..., marker=".", ms=ms, c=mc)
end