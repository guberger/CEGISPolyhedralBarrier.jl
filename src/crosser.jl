struct CrosserProblem
    N::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    afs_inside::Vector{AffForm}
    afs_outside::Vector{AffForm}
    xmax::Float64
    isint::Bool
end

function _optimize_check_crosser!(model, isint)
    optimize!(model)
    isfeasible = primal_status(model) == FEASIBLE_POINT &&
                 (isint || dual_status(model) == FEASIBLE_POINT) &&
                 termination_status(model) == OPTIMAL
    @assert isfeasible || (primal_status(model) == NO_SOLUTION &&
                           termination_status(model) == INFEASIBLE)
    return isfeasible
end

function find_crosser(prob::CrosserProblem, solver)
    model = solver()
    x = @variable(model, [1:prob.N],
                  lower_bound=-prob.xmax,
                  upper_bound=prob.xmax,
                  integer=prob.isint)
    r = @variable(model, upper_bound=10 * prob.xmax)
    y = prob.A * x + prob.b

    for af in prob.afs_inside
        @constraint(model, _eval(af, x) ≤ 0)
    end
    for af in prob.afs_outside
        @constraint(model, _eval(af, y) ≥ r)
    end

    @objective(model, Max, r)

    isfeasible = _optimize_check_crosser!(model, prob.isint)

    xopt = isfeasible ? value.(x) : fill(NaN, prob.N)
    ropt = isfeasible ? value(r) : -Inf

    return xopt, ropt, isfeasible
end