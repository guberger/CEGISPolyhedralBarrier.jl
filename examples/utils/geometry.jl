using LinearAlgebra
using PyCall
const spatial = pyimport_conda("scipy.spatial", "scipy")
const optimize = pyimport_conda("scipy.optimize", "scipy")

function compute_vertices_hrep(halfspaces)
    A = zeros(length(halfspaces), 3)
    A_ub = zeros(length(halfspaces), 3)
    b_ub = zeros(length(halfspaces))
    for (i, h) in enumerate(halfspaces)
        for k = 1:2
            A[i, k] = h[1][k]
            A_ub[i, k] = h[1][k]
        end
        A[i, 3] = h[2]
        A_ub[i, 3] = norm(h[1])
        b_ub[i] = -h[2]
    end
    c = [0, 0, -1]
    bounds = ((nothing, nothing), (nothing, nothing), (nothing, 100))
    res = optimize.linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds)
    @assert res["success"] && res["status"] == 0
    res["fun"] > 0 && return Vector{Float64}[]
    point = res["x"][1:2]
    hs = spatial.HalfspaceIntersection(A, point)
    points = collect.(eachrow(hs.intersections))
    ch = spatial.ConvexHull(points)
    return [ch.points[i + 1, :] for i in ch.vertices]
end

function compute_vertices_vrep(points)
    ch = spatial.ConvexHull(points)
    return [ch.points[i + 1, :] for i in ch.vertices]
end