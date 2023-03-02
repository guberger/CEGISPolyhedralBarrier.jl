function compute_af(
        points_inside::Vector{<:AbstractVector},
        points_image::Vector{<:AbstractVector},
        point_outside::AbstractVector,
        βmax, N, solver
    )
    model = solver()
    a = @variable(model, [1:N], lower_bound=-1, upper_bound=1)
    β = @variable(model, lower_bound=-βmax, upper_bound=βmax)
    r = @variable(model, upper_bound=10)
    af = AffForm(a, β)

    for point in points_inside
        @constraint(model, _eval(af, point) ≤ 0)
    end

    for point in points_image
        @constraint(model, _eval(af, point) + r ≤ 0)
    end

    @constraint(model, _eval(af, point_outside) - r ≥ 0)

    @objective(model, Max, r)

    optimize!(model)

    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT

    return AffForm(value.(af.a), value(af.β)), value(r)
end