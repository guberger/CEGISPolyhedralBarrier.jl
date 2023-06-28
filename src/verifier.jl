struct CexKey
    q::Int # piece
    i::Int # gf_out
end

struct CexVal
    state::State
    r::Float64
end

struct VerifierProblem
    N::Int
    pieces::Vector{Piece}
    gfs_inv::Vector{GenForm}
    gfs_bf::Vector{GenForm}
    gfs_out::Vector{GenForm}
    cexs::Dict{CexKey,CexVal}
    keys_todo::Vector{CexKey}
    afs_inside::Vector{AffForm}
    afs_outside::Vector{AffForm}
end

is_feasible(val, gf) = val.state.loc != gf.loc || _eval(gf.af, val.state.x) ≤ 0

function add_keys_bf_infeasible!(prob::VerifierProblem, indices)
    for (key, val) in prob.cexs
        for i in indices
            if !is_feasible(val, prob.gfs_bf[i])
                push!(prob.keys_todo, key)
                break
            end
        end
    end
end

function add_keys_out_new!(prob::VerifierProblem, indices)
    for (q, piece) in enumerate(prob.pieces)
        for i in indices
            if piece.loc_dst == prob.gfs_out[i].loc
                push!(prob.keys_todo, CexKey(q, i))
            end
        end
    end
end

function make_crossing_problem!(prob, piece, af_out)
    empty!(prob.afs_inside)
    for af in piece.afs_dom
        push!(prob.afs_inside, af)
    end
    for gf in prob.gfs_inv
        if gf.loc == piece.loc_src
            push!(prob.afs_inside, gf.af)
        end
    end
    for gf in prob.gfs_bf
        if gf.loc == piece.loc_src
            push!(prob.afs_inside, gf.af)
        end
    end
    empty!(prob.afs_outside)
    push!(prob.afs_outside, af_out)
    return CrossingProblem(prob.N, piece.A, piece.b,
                           prob.afs_inside, prob.afs_outside)
end

function update_cexs!(prob::VerifierProblem, xmax, solver; int=false)
    for key in prob.keys_todo
        piece, gf_out = prob.pieces[key.q], prob.gfs_out[key.i]
        @assert piece.loc_dst == gf_out.loc
        cross_prob = make_crossing_problem!(prob, piece, gf_out.af)
        x, r, isfeas = find_crosser(cross_prob, xmax, solver, int=int)
        if isfeas
            prob.cexs[key] = CexVal(State(piece.loc_src, x), r)
        else
            delete!(prob.cexs, key)
        end
    end
end