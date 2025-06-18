using LinearAlgebra
using JuMP
using Gurobi
using Plots
using Polyhedra
using CEGISPolyhedralBarrier

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

_isin(afs, x, tol) = all(af -> CPB._eval(af, x) ≤ tol*norm(af.a), afs)

function _next_states!(next_states, pieces, current_state, tol_dom)
    loc, x = current_state.loc, current_state.x
    for piece in pieces
        if piece.loc_src == loc && _isin(piece.afs_dom, x, tol_dom)
            push!(next_states, State(piece.loc_dst, piece.A*x + piece.b))
        end
    end
end

function build_trajectories(pieces, states_init, nstep, tol_dom)
    trajectories = [[state] for state in states_init]
    next_trajectories = Vector{State}[]
    next_states = State[]
    for _ = 1:nstep
        empty!(next_trajectories)
        for traj in trajectories
            empty!(next_states)
            _next_states!(next_states, pieces, traj[end], tol_dom)
            @assert !isempty(next_states)
            for next_state in next_states
                traj_ext = copy(traj)
                push!(traj_ext, next_state)
                push!(next_trajectories, traj_ext)
            end
        end
        trajectories, next_trajectories = next_trajectories, trajectories
    end
    return trajectories
end

function _plot_poly2D!(ax, P; kwargs...)
    plot!(ax, P; kwargs...)
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

function plot_level2D!(ax, afs::Vector{AffForm}, lims; lw=1, kwargs...)
    @assert all(af -> size(af.a) == (2,), afs)
    P = build_polyhedron(afs, lims, 2)
    _plot_poly2D!(ax, P; lw=lw, kwargs...)
end

function plot_level2D!(ax, gfs::Vector{GenForm}, loc::Int, lims;
                       lw=1, kwargs...)
    afs = extract_afs(gfs, loc)
    plot_level2D!(ax, afs, lims; lw=lw, kwargs...)
end

function plot_trajectories2D!(ax,
                              trajectories::Vector{Vector{State}},
                              loc::Int;
                              c=:black, markershape=:circle,
                              lw=0.5, ms=2, kwargs...)
    plot_trajectories2D!(ax, trajectories, (1, 2), (loc,);
                         c=c, markershape=markershape,
                         lw=lw, ms=ms, kwargs...)
end

function plot_trajectories2D!(ax,
                              trajectories::Vector{Vector{State}},
                              vars::Tuple{Int,Int},
                              locs;
                              c=:black, markershape=:circle,
                              lw=0.5, ms=2, kwargs...)
    traj_proj = State[]
    for traj in trajectories
        empty!(traj_proj)
        for state in traj
            if state.loc ∈ locs
                push!(traj_proj, state)
            end
        end
        isempty(traj_proj) && continue
        x1s = [state.x[vars[1]] for state in traj_proj]
        x2s = [state.x[vars[2]] for state in traj_proj]
        plot!(ax, x1s, x2s;
              c=c, markershape=markershape,
              lw=lw, ms=ms, kwargs...)
    end
end

function plot_timeseries!(ax,
                          trajectories::Vector{Vector{State}},
                          var::Int,
                          c=:black, markershape=:circle,
                          lw=0.5, ms=2, kwargs...)
    for traj in trajectories
        xs = [state.x[var] for state in traj]
        ts = collect(0:(length(xs) - 1))
        plot!(ax, ts, xs;
              c=c, markershape=markershape,
              lw=lw, ms=ms, kwargs...)
    end
end