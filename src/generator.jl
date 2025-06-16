struct GeneratorProblem
    N::Int
    gfs::Vector{GenForm}
    gfs_safe::Vector{GenForm}
    indices_new::Vector{Int}
    states_inside::Vector{State}
    states_image::Vector{State}
    edges_unknown::Vector{Edge}
    edges_unknown_new::Vector{Edge}
    xs_inside::Vector{Vector{Float64}}
    xs_inside_margin::Vector{Vector{Float64}}
    xs_outside_margin::Vector{Vector{Float64}}
    ϵ::Float64
    βmax::Float64
    isint::Bool
end

function separation_problem!(prob, state_out)
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
    empty!(prob.xs_outside_margin)
    push!(prob.xs_outside_margin, state_out.x)
    return SeparationProblem(prob.N,
                             prob.xs_inside,
                             prob.xs_inside_margin,
                             prob.xs_outside_margin,
                             prob.βmax, prob.isint)
end

function _inside_margin(prob, gf_safe)
    r_inside::Float64 = Inf
    for state in prob.states_inside
        if state.loc == gf_safe.loc
            r_inside = min(r_inside, -_eval(gf_safe.af, state.x))
        end
    end
    r_image::Float64 = Inf
    for state in prob.states_image
        if state.loc == gf_safe.loc
            r_image = min(r_image, -_eval(gf_safe.af, state.x))
        end
    end
    return r_inside, r_image
end

function add_gfs_safe!(prob)
    for gf_safe in prob.gfs_safe
        r_inside, r_image = _inside_margin(prob, gf_safe)
        if r_inside < 0 || r_image < prob.ϵ
            return false
        end
        r = min(r_inside, r_image - prob.ϵ)
        βnew = min(gf_safe.af.β + r, prob.βmax)
        af = AffForm(gf_safe.af.a, βnew)
        push!(prob.gfs, GenForm(gf_safe.loc, af))
        push!(prob.indices_new, length(prob.gfs))
    end
    return true
end

function update_generator!(prob::GeneratorProblem, solver)
    isreset = false
    empty!(prob.indices_new)
    if isempty(prob.gfs)
        isreset = true
        issuccess = add_gfs_safe!(prob)
        if !issuccess
            return isreset, false
        end
    end
    while !isempty(prob.edges_unknown_new)
        while !isempty(prob.edges_unknown_new)
            edge = pop!(prob.edges_unknown_new)
            sep_prob = separation_problem!(prob, edge.src)
            af, r = find_separator(sep_prob, solver)
            if r < prob.ϵ
                isreset = true
                empty!(prob.gfs)
                empty!(prob.indices_new)
                # push!(prob.states_inside, edge.src) # seems to work well
                push!(prob.states_image, edge.dst)
                append!(prob.edges_unknown_new, prob.edges_unknown)
                empty!(prob.edges_unknown)
                issuccess = add_gfs_safe!(prob)
                if !issuccess
                    return isreset, false
                end
                break # starts with outside first
            else
                push!(prob.gfs, GenForm(edge.src.loc, af))
                push!(prob.indices_new, length(prob.gfs))
                push!(prob.edges_unknown, edge)
            end
        end
    end
    return isreset, true
end