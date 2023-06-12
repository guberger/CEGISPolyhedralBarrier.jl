struct SeparationProblem
    N::Int
    grid_inside::Grid
    grid_inside_margin::Grid
    grid_outside::Grid
    grid_outside_margin::Grid
end

struct AffFormVar
    a::Vector{VariableRef}
    β::VariableRef
end
_eval(af::AffFormVar, point) = dot(af.a, point) + af.β

function _optimize_check_separator!(model)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT
    @assert dual_status(model) == FEASIBLE_POINT
    return nothing
end

function find_separator(prob::SeparationProblem, βmax, solver)
    model = solver()
    a = @variable(model, [1:prob.N], lower_bound=-1, upper_bound=1)
    β = @variable(model, lower_bound=-βmax, upper_bound=βmax)
    r = @variable(model, upper_bound=10)
    af = AffFormVar(a, β)

    for point in prob.grid_inside.points
        @constraint(model, _eval(af, point) ≤ 0)
    end
    for point in prob.grid_inside_margin.points
        @constraint(model, _eval(af, point) + r ≤ 0)
    end
    for point in prob.grid_outside.points
        @constraint(model, _eval(af, point) ≥ 0)
    end
    for point in prob.grid_outside_margin.points
        @constraint(model, _eval(af, point) - r ≥ 0)
    end

    @objective(model, Max, r)

    _optimize_check_separator!(model)

    return AffForm(value.(af.a), value(af.β)), value(r)
end