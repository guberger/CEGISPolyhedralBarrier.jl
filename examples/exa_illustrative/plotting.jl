function _to_vector_ineq(poly::Polyhedron, nvar)
    A = zeros(length(poly.halfspaces), nvar)
    b = zeros(length(poly.halfspaces))
    for (i, h) in enumerate(poly.halfspaces)
        for k = 1:nvar
            A[i, k] = h.a[k]
        end
        b[i] = -h.β
    end
    return A, b
end

function _to_vector_ineq(pf::PolyFunc, nvar)
    A = zeros(length(pf.afs), nvar)
    b = zeros(length(pf.afs))
    for (i, af) in enumerate(pf.afs)
        for k = 1:nvar
            A[i, k] = af.lin[k]
        end
        b[i] = -af.off
    end
    return A, b
end

function plot_hrep!(
        ax, A, b; fc="blue", fa=0.5, ec="blue", ew=2.0
    )
    verts = compute_vertices_hrep(A, b)
    isempty(verts) && return
    polylist = matplotlib.collections.PolyCollection((verts,))
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function plot_vrep!(
        ax, points; fc="blue", fa=0.5, ec="blue", ew=2.0
    )
    isempty(points) && return
    verts = compute_vertices_vrep(points)
    polylist = matplotlib.collections.PolyCollection((verts,))
    fca = matplotlib.colors.colorConverter.to_rgba(fc, alpha=fa)
    polylist.set_facecolor(fca)
    polylist.set_edgecolor(ec)
    polylist.set_linewidth(ew)
    ax.add_collection(polylist)
end

function plot_point!(ax, point; mc="blue", ms=15)
    ax.plot(point..., marker=".", ms=ms, c=mc)
end