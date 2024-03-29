struct CrossingProblem
    N::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    afs_inside::Vector{AffForm}
    afs_outside::Vector{AffForm}
end

function _optimize_check_counterexample!(model, int)
    optimize!(model)
    isfeasible = primal_status(model) == FEASIBLE_POINT &&
                 (int || dual_status(model) == FEASIBLE_POINT) &&
                 termination_status(model) == OPTIMAL
    @assert isfeasible || (primal_status(model) == NO_SOLUTION &&
                           termination_status(model) == INFEASIBLE)
    return isfeasible
end

function find_crosser(prob::CrossingProblem, xmax, solver; int=false)
    model = solver()
    x = @variable(model, [1:prob.N],
                  lower_bound=-xmax, upper_bound=xmax, integer=int)
    r = @variable(model, upper_bound=10*xmax)
    y = prob.A*x + prob.b

    for af in prob.afs_inside
        @constraint(model, _eval(af, x) ≤ 0)
    end
    for af in prob.afs_outside
        @constraint(model, _eval(af, y) ≥ norm(af.a, Inf)*r)
    end

    @objective(model, Max, r)

    isfeasible = _optimize_check_counterexample!(model, int)

    xopt = isfeasible ? value.(x) : fill(NaN, prob.N)
    ropt = isfeasible ? value(r) : -Inf

    return xopt, ropt, isfeasible
end