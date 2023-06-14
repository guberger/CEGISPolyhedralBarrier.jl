struct CrossingProblem
    N::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    afs_inside_fragile::Vector{AffForm}
    afs_inside_margin::Vector{AffForm}
    afs_outside::Vector{AffForm}
    λ::Float64
end

function _optimize_check_counterexample!(model)
    optimize!(model)
    isfeasible = primal_status(model) == FEASIBLE_POINT &&
                 dual_status(model) == FEASIBLE_POINT &&
                 termination_status(model) == OPTIMAL
    @assert isfeasible || (primal_status(model) == NO_SOLUTION &&
                           termination_status(model) == INFEASIBLE)
    return isfeasible
end

function find_crosser(prob::CrossingProblem, xmax, solver)
    model = solver()
    x = @variable(model, [1:prob.N], lower_bound=-xmax, upper_bound=xmax)
    r = @variable(model, upper_bound=10*xmax)
    y = prob.A*x + prob.b

    for af in prob.afs_inside_fragile
        @constraint(model, _eval(af, x) ≤ 0)
    end
    for af in prob.afs_inside_margin
        @constraint(model, _eval(af, x) ≤ -margin(af, 0, r))
    end
    for af in prob.afs_outside
        @constraint(model, _eval(af, y) ≥ prob.λ*margin(af, 0, r))
    end

    @objective(model, Max, r)

    isfeasible = _optimize_check_counterexample!(model)

    xopt = isfeasible ? value.(x) : fill(NaN, prob.N)
    ropt = isfeasible ? value(r) : -Inf

    return xopt, ropt, isfeasible
end