struct CexKey
    q::Int # piece
    i::Int # gf_bf
end

struct CexVal
    state::State
    r::Float64
end

struct VerifierProblem
    N::Int
    pieces::Vector{Piece}
    gfs_bf::Vector{GenForm}
    gfs_safe::Vector{GenForm}
    gfs_inv::Vector{GenForm}
    cexs::Dict{CexKey,CexVal}
    keys_todo::Vector{CexKey}
    afs_inside::Vector{AffForm}
    afs_inside_margin::Vector{AffForm}
    afs_outside::Vector{AffForm}
    afs_outside_margin::Vector{AffForm}
end

function add_infeasible_keys!(prob::VerifierProblem, gfs)
    for (key, val) in prob.cexs
        state = val.state
        for gf in gfs
            if state.loc == gf.loc && _eval(gf.af, state.x) > 0
                push!(prob.keys_todo, key)
                break
            end
        end
    end
end

function add_gfs_keys!(prob::VerifierProblem, gfs, indices)
    for (q, piece) in enumerate(prob.pieces)
        for i in indices
            if piece.loc_dst == gfs[i].loc
                push!(prob.keys_todo, CexKey(q, i))
            end
        end
    end
end

function prepare_crossing_problem!(prob, piece)
    empty!(prob.afs_inside)
    for gf in prob.gfs_inv
        if gf.loc == piece.loc_src
            push!(prob.afs_inside, gf.af)
        end
    end
    for af in piece.afs_dom
        push!(prob.afs_inside, af)
    end
    empty!(prob.afs_inside_margin)
    for gf in prob.gfs_safe
        if gf.loc == piece.loc_src
            push!(prob.afs_inside_margin, gf.af)
        end
    end
    for gf in prob.gfs_bf
        if gf.loc == piece.loc_src
            push!(prob.afs_inside_margin, gf.af)
        end
    end
    empty!(prob.afs_outside)
    empty!(prob.afs_outside_margin)
    return CrossingProblem(prob.N, piece.A, piece.b,
                           prob.afs_inside, prob.afs_inside_margin,
                           prob.afs_outside, prob.afs_outside_margin)
end

function update_cexs_safe!(prob::VerifierProblem, xmax, solver)
    for key in prob.keys_todo
        piece, gf_out = prob.pieces[key.q], prob.gfs_safe[key.i]
        @assert piece.loc_dst == gf_out.loc
        cross_prob = prepare_crossing_problem!(prob, piece)
        push!(cross_prob.afs_outside, gf_out.af)
        x, r, isfeasible = find_crosser(cross_prob, xmax, solver)
        if isfeasible
            prob.cexs[key] = CexVal(State(piece.loc_src, x), r)
        end
    end
end

function update_cexs_cont!(prob::VerifierProblem, xmax, solver)
    for key in prob.keys_todo
        piece, gf_out = prob.pieces[key.q], prob.gfs_bf[key.i]
        @assert piece.loc_dst == gf_out.loc
        cross_prob = prepare_crossing_problem!(prob, piece)
        push!(cross_prob.afs_outside_margin, gf_out.af)
        x, r, isfeasible = find_crosser(cross_prob, xmax, solver)
        if isfeasible
            prob.cexs[key] = CexVal(State(piece.loc_src, x), r)
        end
    end
end