struct CrossingProblem
    N::Int
    afs_inside::Vector{AffForm}
    afs_inside_margin::Vector{AffForm}
    afs_outside::Vector{AffForm}
    afs_outside_margin::Vector{AffForm}
end
empty_cross_problem(N) = CrossingProblem(N, ntuple(i -> AffForm[], Val(4))...)

function _optimize_check_counterexample!(model)
    optimize!(model)
    flag = primal_status(model) == FEASIBLE_POINT &&
           dual_status(model) == FEASIBLE_POINT &&
           termination_status(model) == OPTIMAL
    @assert flag || (primal_status(model) == NO_SOLUTION &&
                     termination_status(model) == INFEASIBLE)
    return flag
end

function find_crosser(prob::CrossingProblem, A, b, xmax, solver)
    model = solver()
    x = @variable(model, [1:prob.N], lower_bound=-xmax, upper_bound=xmax)
    r = @variable(model, upper_bound=10*xmax)

    for af in prob.afs_inside
        @constraint(model, _eval(af, x) ≤ 0)
    end
    for af in prob.afs_inside_margin
        @constraint(model, _eval(af, x) + margin(af, 0, r) ≤ 0)
    end
    for af in prob.afs_outside
        @constraint(model, _eval(af, A*x + b) ≥ 0)
    end
    for af in prob.afs_outside_margin
        @constraint(model, _eval(af, A*x + b) ≥ margin(af, 0, r))
    end

    @objective(model, Max, r)

    flag = _optimize_check_counterexample!(model)

    xopt = flag ? value.(x) : fill(NaN, N)
    ropt = flag ? value(r) : -Inf

    return xopt, ropt, flag
end

function _reset_crossing_problem!(prob)
    empty!(prob.afs_inside)
    empty!(prob.afs_inside_margin)
    empty!(prob.afs_outside)
    empty!(prob.afs_outside_margin)
end

struct VerifierSafeProblem
    N::Int
    sys::System
    mpf_safe::MultiPolyFunc
    mpf_inv::MultiPolyFunc
    mpf_BF::MultiPolyFunc
end

function find_counterexample(prob::VerifierSafeProblem, xmax, solver)
    cross_prob = empty_cross_problem(prob.N)
    xopt::Vector{Float64} = fill(NaN, prob.N)
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in prob.sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        _reset_crossing_problem!(cross_prob)
        append!(cross_prob.afs_inside, pf_dom.afs)
        append!(cross_prob.afs_inside, prob.mpf_inv.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_safe.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_BF.pfs[piece.loc1].afs)
        for af in prob.mpf_safe.pfs[piece.loc2].afs
            empty!(cross_prob.afs_outside)
            push!(cross_prob.afs_outside, af)
            x, r, flag = find_crosser(cross_prob, A, b, xmax, solver
            )
            if flag && r > ropt
                xopt = x
                ropt = r
                locopt = piece.loc1
            end
        end
    end
    return xopt, ropt, locopt
end

struct VerifierBFProblem
    N::Int
    sys::System
    mpf_safe::MultiPolyFunc
    mpf_inv::MultiPolyFunc
    mpf_BF::MultiPolyFunc
end

function find_counterexample(prob::VerifierBFProblem, xmax, solver)
    cross_prob = empty_cross_problem(prob.N)
    xopt::Vector{Float64} = fill(NaN, prob.N)
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in prob.sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        _reset_crossing_problem!(cross_prob)
        append!(cross_prob.afs_inside, pf_dom.afs)
        append!(cross_prob.afs_inside, prob.mpf_inv.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_safe.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_BF.pfs[piece.loc1].afs)
        for af in prob.mpf_BF.pfs[piece.loc2].afs
            empty!(cross_prob.afs_outside_margin)
            push!(cross_prob.afs_outside_margin, af)
            x, r, flag = find_crosser(cross_prob, A, b, xmax, solver
            )
            if flag && r > ropt
                xopt = x
                ropt = r
                locopt = piece.loc1
            end
        end
    end
    return xopt, ropt, locopt
end