ContPiece = CEGISPolyhedralBarrier.PieceCont
PolyFunc = CEGISPolyhedralBarrier.PolyFunc
MultiPolyFunc = CEGISPolyhedralBarrier.MultiPolyFunc
Polyhedron = CEGISPolyhedralBarrier.Polyhedron
PosEvidence = CEGISPolyhedralBarrier.PosEvidence
LieContEvidence = CEGISPolyhedralBarrier.LieContEvidence
_eval = CEGISPolyhedralBarrier._eval
_norm(pf::PolyFunc, x) = maximum(lf -> _eval(lf, x), pf.lfs)

function plot_field!(ax, piece::ContPiece, xlims, ylims, ngrid; c="gray")
    x1_grid = range(xlims..., length=ngrid)
    x2_grid = range(ylims..., length=ngrid)
    X = map(x -> [x...], Iterators.product(x1_grid, x2_grid))
    X1 = getindex.(X, 1)
    X2 = getindex.(X, 2)
    D = Matrix{Vector{Float64}}(undef, ngrid, ngrid)
    for (k, x) in enumerate(X)
        if x ∈ piece.domain
            D[k] = piece.A*x
        else
            D[k] = [NaN, NaN]
        end
    end
    D1 = getindex.(D, 1)
    D2 = getindex.(D, 2)
    ax.quiver(X1, X2, D1, D2, color=c)
end

function plot_level!(
        ax, pf::PolyFunc, radmax;
        fc="gold", fa=0.5, ec="gold", ew=2.0
    )
    p = Polyhedron()
    for lf in pf.lfs
        CPB.add_halfspace!(p, lf.lin, -1)
    end
    verts = compute_vertices_2d(p, zeros(2))
    verts_radius = maximum(vert -> norm(vert, Inf), verts)
    scaling = radmax/verts_radius
    verts_scaled = map(vert -> vert*scaling, verts)
    polylist = matplotlib.collections.PolyCollection([verts_scaled])
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
    return scaling
end

function plot_evids!(
        ax, evids::Vector{PosEvidence}, mpf::MultiPolyFunc, level;
        mc="blue", ms=15
    )
    for evid in evids
        point_norm = _norm(mpf.pfs[evid.loc], evid.point)
        point_scaled = evid.point*level/point_norm
        ax.plot(point_scaled..., marker=".", ms=ms, c=mc)
    end
end

_get_deriv(evid::LieContEvidence) = (evid.point2 - evid.point1)/evid.τ

function plot_evids!(
        ax, evids::Vector{LieContEvidence},
        mpf::MultiPolyFunc, level, deriv_length;
        mc="blue", ms=15, lc="green", lw=2.5
    )
    derivs_norm = -Inf
    for evid in evids
        point1_norm = _norm(mpf.pfs[evid.loc], evid.point1)
        deriv = _get_deriv(evid)
        derivs_norm = max(derivs_norm, norm(deriv)/point1_norm)
    end
    for evid in evids
        point1_norm = _norm(mpf.pfs[evid.loc], evid.point1)
        point1_scaled = evid.point1*level/point1_norm
        ax.plot(point1_scaled..., marker=".", ms=ms, c=mc)
        deriv = _get_deriv(evid)
        deriv_scaled = deriv/(point1_norm*derivs_norm)
        point2_scaled = point1_scaled + deriv_length*deriv_scaled
        ax.plot(
            (point1_scaled[1], point2_scaled[1]),
            (point1_scaled[2], point2_scaled[2]), c=lc, lw=lw
        )
    end
end

function plot_traj!(ax, sys, x0, dt, nstep; c="purple", ms=15, lw=2.5)
    ax.plot(x0..., marker=".", ms=ms, c=c)
    x_seq = Vector{Vector{Float64}}(undef, nstep)
    x_seq[1] = x0
    for i = 2:nstep
        x = x_seq[i - 1]
        next = false
        for piece in sys.disc_pieces
            next && break
            if x ∈ piece.domain
                x = piece.A*x
                next = true
            end
        end
        for piece in sys.cont_pieces
            next && break
            if x ∈ piece.domain
                x = exp(piece.A*dt)*x
                next = true
            end
        end
        x_seq[i] = x
    end
    ax.plot(getindex.(x_seq, 1), getindex.(x_seq, 2), lw=lw, c=c)
end