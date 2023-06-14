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
    gfs_safe::Vector{GenForm}
    gfs_bf::Vector{GenForm}
    gfs_out::Vector{GenForm}
    cexs::Dict{CexKey,CexVal}
    keys_todo::Vector{CexKey}
    afs_inside_fragile::Vector{AffForm}
    afs_inside_margin::Vector{AffForm}
    afs_outside::Vector{AffForm}
    λ::Float64
end

_is_feasible(val::CexVal, gf::GenForm) =
    val.state.loc != gf.loc || _eval(gf.af, val.state.x) ≤ -margin(gf.af, 0, val.r)

function add_keys_bf_infeasible!(prob::VerifierProblem, indices)
    for (key, val) in prob.cexs
        for i in indices
            if !_is_feasible(val, prob.gfs_bf[i])
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
    empty!(prob.afs_inside_fragile)
    for af in piece.afs_dom
        push!(prob.afs_inside_fragile, af)
    end
    for gf in prob.gfs_inv
        if gf.loc == piece.loc_src
            push!(prob.afs_inside_fragile, gf.af)
        end
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
    push!(prob.afs_outside, af_out)
    return CrossingProblem(prob.N, piece.A, piece.b,
                           prob.afs_inside_fragile,
                           prob.afs_inside_margin,
                           prob.afs_outside, prob.λ)
end

function update_cexs!(prob::VerifierProblem, xmax, solver)
    for key in prob.keys_todo
        piece, gf_out = prob.pieces[key.q], prob.gfs_out[key.i]
        @assert piece.loc_dst == gf_out.loc
        cross_prob = make_crossing_problem!(prob, piece, gf_out.af)
        x, r, isfeasible = find_crosser(cross_prob, xmax, solver)
        if isfeasible
            prob.cexs[key] = CexVal(State(piece.loc_src, x), r)
        end
    end
end