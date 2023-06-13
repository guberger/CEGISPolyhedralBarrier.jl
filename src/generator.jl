struct GeneratorProblem
    N::Int
    mpf::Vector{PolyFunc}
    mgrid_inside::Vector{Grid}
    mgrid_image::Vector{Grid}
    graph_unknown::Graph
    graph_unknown_new::Graph
    graph_outside::Graph
    graph_outside_new::Graph
    graph_temp::Graph
    ϵ::Float64
end

function update_generator!(prob::GeneratorProblem, βmax, solver)
    while !isempty(prob.graph_outside_new.links) ||
          !isempty(prob.graph_unknown_new.links)
        while !isempty(prob.graph_outside_new.links)
            # print(".")
            link = pop!(prob.graph_outside_new.links)
            sep_prob = SeparationProblem(prob.N,
                                         prob.mgrid_inside[link.loc_pre],
                                         prob.mgrid_image[link.loc_pre],
                                         empty_grid(),
                                         Grid([link.point_pre]))
            af, r = find_separator(sep_prob, βmax, solver)
            if r < prob.ϵ
                return false
            else
                push!(prob.mpf[link.loc_pre].afs, af)
                push!(prob.graph_outside.links, link)
            end
        end
        while !isempty(prob.graph_unknown_new.links)
            # print("+")
            link = pop!(prob.graph_unknown_new.links)
            sep_prob = SeparationProblem(prob.N,
                                         prob.mgrid_inside[link.loc_pre],
                                         prob.mgrid_image[link.loc_pre],
                                         empty_grid(),
                                         Grid([link.point_pre]))
            af, r = find_separator(sep_prob, βmax, solver)
            if r < prob.ϵ
                push!(prob.mgrid_inside[link.loc_pre].points, link.point_pre)
                push!(prob.mgrid_image[link.loc_post].points, link.point_post)
                empty!(prob.mpf[link.loc_pre].afs)
                empty!(prob.mpf[link.loc_post].afs)
                empty!(prob.graph_temp.links)
                while !isempty(prob.graph_outside.links)
                    link2 = pop!(prob.graph_outside.links)
                    if link2.loc_pre == link.loc_pre ||
                        link2.loc_pre == link.loc_post
                        push!(prob.graph_outside_new.links, link2)
                    else
                        push!(prob.graph_temp.links, link2)
                    end
                end
                for link2 in prob.graph_temp.links
                    push!(prob.graph_outside.links, link2)
                end
                empty!(prob.graph_temp.links)
                while !isempty(prob.graph_unknown.links)
                    link2 = pop!(prob.graph_unknown.links)
                    if link2.loc_pre == link.loc_pre ||
                        link2.loc_pre == link.loc_post
                        push!(prob.graph_unknown_new.links, link2)
                    else
                        push!(prob.graph_temp.links, link2)
                    end
                end
                for link2 in prob.graph_temp.links
                    push!(prob.graph_unknown.links, link2)
                end
                break # starts with outside first
            else
                push!(prob.mpf[link.loc_pre].afs, af)
                push!(prob.graph_unknown.links, link)
            end
        end
    end
    return true
end