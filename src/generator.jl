struct GeneratorProblem
    N::Int
    gfs::Vector{GenForm}
    indices_new::Vector{Int}
    states_inside::Vector{State}
    states_image::Vector{State}
    links_unknown::Vector{Link}
    links_unknown_new::Vector{Link}
    states_outside::Vector{State}
    states_outside_new::Vector{State}
    xs_inside::Vector{Vector{Float64}}
    xs_inside_margin::Vector{Vector{Float64}}
    xs_outside::Vector{Vector{Float64}}
    xs_outside_margin::Vector{Vector{Float64}}
    ϵ::Float64
end

function make_separation_problem!(prob, state_out)
    empty!(prob.xs_inside)
    for state in prob.states_inside
        if state.loc == state_out.loc
            push!(prob.xs_inside, state.x)
        end
    end
    empty!(prob.xs_inside_margin)
    for state in prob.states_image
        if state.loc == state_out.loc
            push!(prob.xs_inside_margin, state.x)
        end
    end
    empty!(prob.xs_outside)
    empty!(prob.xs_outside_margin)
    push!(prob.xs_outside_margin, state_out.x)
    return SeparationProblem(prob.N,
                             prob.xs_inside, prob.xs_inside_margin,
                             prob.xs_outside, prob.xs_outside_margin)
end

function update_generator!(prob::GeneratorProblem, βmax, solver)
    isreset = false
    empty!(prob.indices_new)
    while !isempty(prob.states_outside_new) || !isempty(prob.links_unknown_new)
        while !isempty(prob.states_outside_new)
            state = pop!(prob.states_outside_new)
            sep_prob = make_separation_problem!(prob, state)
            af, r = find_separator(sep_prob, βmax, solver)
            if r < prob.ϵ
                push!(prob.states_outside_new, state)
                return isreset, false
            else
                push!(prob.gfs, GenForm(state.loc, af))
                push!(prob.indices_new, length(prob.gfs))
                push!(prob.states_outside, state)
            end
        end
        while !isempty(prob.links_unknown_new)
            link = pop!(prob.links_unknown_new)
            sep_prob = make_separation_problem!(prob, link.src)
            af, r = find_separator(sep_prob, βmax, solver)
            if r < prob.ϵ
                isreset = true
                empty!(prob.gfs)
                empty!(prob.indices_new)
                push!(prob.states_inside, link.src)
                push!(prob.states_image, link.dst)
                append!(prob.links_unknown_new, prob.links_unknown)
                empty!(prob.links_unknown)
                append!(prob.states_outside_new, prob.states_outside)
                empty!(prob.states_outside)
                break # starts with outside first
            else
                push!(prob.gfs, GenForm(link.src.loc, af))
                push!(prob.indices_new, length(prob.gfs))
                push!(prob.links_unknown, link)
            end
        end
    end
    return isreset, true
end