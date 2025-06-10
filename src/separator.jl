struct SeparationProblem
    N::Int
    xs_inside::Vector{Vector{Float64}}
    xs_inside_margin::Vector{Vector{Float64}}
    xs_outside::Vector{Vector{Float64}}
    ϵ::Float64
    βmax::Float64
    isint::Bool
end

struct AffFormVar
    a::Vector{VariableRef}
    β::VariableRef
end
_eval(af::AffFormVar, x) = dot(af.a, x) + af.β

function _optimize_check_separator!(model, isint)
    optimize!(model)
    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT
    @assert isint || dual_status(model) == FEASIBLE_POINT
    return nothing
end

function find_separator(prob::SeparationProblem, solver)
    model = solver()
    a = @variable(model, [1:prob.N],
                  lower_bound=-1, upper_bound=1, integer=prob.isint)
    β = @variable(model, lower_bound=-prob.βmax, upper_bound=prob.βmax)
    r = @variable(model, upper_bound=10) # 10 is an arbitrary parameter
    af = AffFormVar(a, β)

    for x in prob.xs_inside
        @constraint(model, _eval(af, x) + r ≤ 0)
    end
    for x in prob.xs_inside_margin
        @constraint(model, _eval(af, x) + r + prob.ϵ ≤ 0)
    end
    for x in prob.xs_outside
        @constraint(model, _eval(af, x) - r ≥ 0)
    end

    @objective(model, Max, r)

    _optimize_check_separator!(model, prob.isint)

    return AffForm(value.(af.a), value(af.β)), value(r)
end