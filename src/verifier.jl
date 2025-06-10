struct Witness
    q::Int # piece
    i::Int # gf_out
    state::State
    r::Float64
    isfeas::Bool
end

struct VerifierProblem
    N::Int
    pieces::Vector{Piece}
    gfs_inv::Vector{GenForm}
    gfs_bf::Vector{GenForm}
    gfs_out::Vector{GenForm}
    witnesses::Vector{Witness}
    keys_todo::Vector{Tuple{Int,Int,Int}}
    afs_inside::Vector{AffForm}
    afs_outside::Vector{AffForm}
    xmax::Float64
    isint::Bool
end

function add_keys_bf_infeasible!(prob::VerifierProblem, indices)
    for (k, witness) in enumerate(prob.witnesses)
        state = witness.state
        for i in indices
            gf = prob.gfs_bf[i]
            if state.loc == gf.loc && _eval(gf.af, state.x) > 0
                push!(prob.keys_todo, (witness.q, witness.i, k))
                break
            end
        end
    end
end

function add_keys_out_new!(prob::VerifierProblem, indices)
    for (q, piece) in enumerate(prob.pieces)
        for i in indices
            if piece.loc_dst == prob.gfs_out[i].loc
                push!(prob.keys_todo, (q, i, -1))
            end
        end
    end
end

function crosser_problem!(prob, piece, af_out)
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
    return CrosserProblem(prob.N, piece.A, piece.b,
                           prob.afs_inside, prob.afs_outside,
                           prob.xmax, prob.isint)
end

function update_verifier!(prob::VerifierProblem, solver)
    for (q, i, k) in prob.keys_todo
        piece, gf_out = prob.pieces[q], prob.gfs_out[i]
        @assert piece.loc_dst == gf_out.loc
        cross_prob = crosser_problem!(prob, piece, gf_out.af)
        x, r, isfeas = find_crosser(cross_prob, solver)
        witness = Witness(q, i, State(piece.loc_src, x), r, isfeas)
        if k > 0
            prob.witnesses[k] = witness
        else
            push!(prob.witnesses, witness)
        end
    end
end