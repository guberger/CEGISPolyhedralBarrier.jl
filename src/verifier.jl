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
    gfs_new::Vector{GenForm}
    gfs_inside::Vector{GenForm}
    gfs_outside::Vector{GenForm}
    indices_inside_new::Vector{Int}
    indices_outside_new::Vector{Int}
    witnesses::Vector{Witness}
    keys_todo::Vector{Tuple{Int,Int,Int}}
    afs_inside::Vector{AffForm}
    afs_outside::Vector{AffForm}
    xmax::Float64
    isint::Bool
end

function _update_keys!(prob::VerifierProblem)
    # add keys that were inside before, but no longer inside
    for (k, witness) in enumerate(prob.witnesses)
        state = witness.state
        for i in prob.indices_inside_new
            gf = prob.gfs_inside[i]
            if state.loc == gf.loc && _eval(gf.af, state.x) > 0
                push!(prob.keys_todo, (witness.q, witness.i, k))
                break
            end
        end
    end
    # add keys for new outside
    for (q, piece) in enumerate(prob.pieces)
        for i in prob.indices_outside_new
            if piece.loc_dst == prob.gfs_outside[i].loc
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
    for gf in prob.gfs_inside
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

function update_verifier!(prob::VerifierProblem, isreset, solver)
    empty!(prob.keys_todo)
    empty!(prob.indices_inside_new)
    empty!(prob.indices_outside_new)
    if isreset
        empty!(prob.witnesses)
        empty!(prob.gfs_inside)
        empty!(prob.gfs_outside)
        copy!(prob.gfs_inside, prob.gfs_inv)
    end
    for gf in prob.gfs_new
        push!(prob.gfs_inside, gf)
        push!(prob.gfs_outside, gf)
        push!(prob.indices_inside_new, length(prob.gfs_inside))
        push!(prob.indices_outside_new, length(prob.gfs_outside))
    end
    _update_keys!(prob)
    for (q, i, k) in prob.keys_todo
        piece, gf_out = prob.pieces[q], prob.gfs_outside[i]
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
    kopt::Int = -1
    rmax::Float64 = -Inf
    for (k, witness) in enumerate(prob.witnesses)
        if witness.isfeas && witness.r > rmax
            rmax = witness.r
            kopt = k
        end
    end
    return kopt, rmax
end