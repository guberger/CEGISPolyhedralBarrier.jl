struct CrossingProblem
    N::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    pf_inside::PolyFunc
    pf_inside_margin::PolyFunc
    pf_outside::PolyFunc
    pf_outside_margin::PolyFunc
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

function find_crosser(prob::CrossingProblem, xmax, solver)
    model = solver()
    x = @variable(model, [1:prob.N], lower_bound=-xmax, upper_bound=xmax)
    r = @variable(model, upper_bound=10*xmax)
    y = prob.A*x + prob.b

    for af in prob.pf_inside.afs
        @constraint(model, _eval(af, x) ≤ 0)
    end
    for af in prob.pf_inside_margin.afs
        @constraint(model, _eval(af, x) ≤ -margin(af, 0, r))
    end
    for af in prob.pf_outside.afs
        @constraint(model, _eval(af, y) ≥ 0)
    end
    for af in prob.pf_outside_margin.afs
        @constraint(model, _eval(af, y) ≥ margin(af, 0, r))
    end

    @objective(model, Max, r)

    flag = _optimize_check_counterexample!(model)

    xopt = flag ? value.(x) : fill(NaN, prob.N)
    ropt = flag ? value(r) : -Inf

    return xopt, ropt, flag
end