struct Separator{N}
    soft_evids::Vector{Point{N}}
    hard_evids::Vector{Point{N}}
    ϵ::Float64
    βmax::Float64
end

Separator{N}(ϵ, βmax) where N = Separator(Point{N}[], Point{N}[], ϵ, βmax)

add_soft_evid!(sep::Separator, point) = push!(sep.soft_evids, point)
add_hard_evid!(sep::Separator, point) = push!(sep.hard_evids, point)

struct _AF{N}
    a::SVector{N,VariableRef}
    β::VariableRef
end
_eval(af::_AF, point) = dot(point, af.a) + af.β

function compute_af(sep::Separator{N}, point::Point{N}, solver) where N
    model = solver()
    a = SVector(ntuple(
        k -> @variable(model, lower_bound=-1, upper_bound=1), Val(N)
    ))
    β = @variable(model, lower_bound=-sep.βmax, upper_bound=sep.βmax)
    r = @variable(model, upper_bound=10)
    af = _AF(a, β)

    for point in sep.soft_evids
        @constraint(model, _eval(af, point) + r ≤ 0)
    end

    for point in sep.hard_evids
        @constraint(model, _eval(af, point) + r + sep.ϵ ≤ 0)
    end

    @constraint(model, _eval(af, point) - r - sep.ϵ ≥ 0)

    @objective(model, Max, r)

    optimize!(model)

    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT

    return AffForm(value.(af.a), value(af.β)), value(r)
end