## Learner

@enum StatusCode begin
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    MAX_ITER_REACHED = 3
end

struct BarrierProblem
    N::Int
    pieces::Vector{Piece}
    gfs_inv::Vector{GenForm}
    gfs_safe::Vector{GenForm}
    states_init::Vector{State}
    ϵ::Float64
    δ::Float64
end

# :tol_dom=1e-8,
# :βmax=1e3,
# :xmax=1e3

function build_generator(prob)
    return GeneratorProblem(
        prob.N,
        GenForm[], # gfs
        Int[], # indices_new
        copy(prob.states_init), # states_inside
        State[], # states_image
        Link[], # links_unknown
        Link[], # links_unknown_new
        State[], # states_outside
        State[], # states_outside_new
        [Vector{Float64}[] for i = 1:4]..., # all xs_...
        prob.ϵ
    )
end

function build_verifier_safe(prob)
    return VerifierProblem(
        prob.N,
        prob.pieces,
        copy(prob.gfs_inv),
        copy(prob.gfs_safe),
        GenForm[], # gfs_bf
        copy(prob.gfs_safe), # gfs_out
        Dict{CexKey,CexVal}(), # cexs
        CexKey[], # keys_todo
        [AffForm[] for i = 1:3]..., # all afs_...
        0.0
    )
end

function build_verifier_cont(prob)
    return VerifierProblem(
        prob.N,
        prob.pieces,
        copy(prob.gfs_inv),
        copy(prob.gfs_safe),
        GenForm[], # gfs_bf
        GenForm[], # gfs_out
        Dict{CexKey,CexVal}(), # cexs
        CexKey[], # keys_todo
        [AffForm[] for i = 1:3]..., # all afs_...
        1.0
    )
end

function reset_verifier!(prob::VerifierProblem)
    empty!(prob.cexs)
    empty!(prob.keys_todo)
    add_keys_out_new!(prob, eachindex(prob.gfs_out))
end

function find_cex_max(cexs::Dict{CexKey,CexVal})
    keyopt::CexKey = CexKey(0, 0)
    rmax::Float64 = -Inf
    for (key, val) in cexs
        if val.r > rmax
            rmax = val.r
            keyopt = key
        end
    end
    return keyopt, rmax
end

## Learn Barrier
function find_barrier(prob::BarrierProblem,
                      iter_max, solver; # LP solver
                      βmax=1e3, xmax=1e3,
                      do_print=true)
    
    gen_prob = build_generator(prob)
    verif_safe_prob = build_verifier_safe(prob)
    reset_verifier!(verif_safe_prob)
    verif_cont_prob = build_verifier_cont(prob)
    reset_verifier!(verif_cont_prob)

    iter = 0
    isfound = false
    issuccess = true
    
    while iter < iter_max && !isfound && issuccess
        iter += 1
        do_print && print("Iter: ", iter)

        # Generation part
        do_print && print("|--- Generate... ")
        isreset, issuccess = update_generator!(gen_prob, βmax, solver)

        if !issuccess
            println("Failed")
            break
        end

        # Verification safe part
        copy!(verif_safe_prob.gfs_bf, gen_prob.gfs)
        copy!(verif_cont_prob.gfs_bf, gen_prob.gfs)
        copy!(verif_cont_prob.gfs_out, gen_prob.gfs)
        if isreset
            do_print && println("Reset")
            reset_verifier!(verif_safe_prob)
            reset_verifier!(verif_cont_prob)
        else
            do_print && println("Expanded")
            add_keys_bf_infeasible!(verif_safe_prob, gen_prob.indices_new)
            add_keys_bf_infeasible!(verif_cont_prob, gen_prob.indices_new)
            add_keys_out_new!(verif_cont_prob, gen_prob.indices_new)
        end

        do_print && print("|--- Verify safe... ")
        update_cexs!(verif_safe_prob, xmax, solver)
        key, r = find_cex_max(verif_safe_prob.cexs)
        if r > -prob.δ
            state = verif_safe_prob.cexs[key].state
            do_print && println("CE found: ", state.loc, ", ", state.x, ", ", r)
            push!(gen_prob.states_outside_new, state)
            continue
        else
            do_print && println("No CE found: ", r)
        end

        do_print && print("|--- Verify cont... ")
        update_cexs!(verif_cont_prob, xmax, solver)
        key, r = find_cex_max(verif_cont_prob.cexs)
        if r > -prob.δ
            state = verif_cont_prob.cexs[key].state
            piece = verif_cont_prob.pieces[key.q]
            do_print && println("CE found: ", state.loc, ", ", state.x, ", ", r)
            state_post = State(piece.loc_dst, piece.A*state.x + piece.b)
            link = Link(state, state_post)
            push!(gen_prob.links_unknown_new, link)
            continue
        else
            do_print && println("No CE found: ", r)
        end

        # Passed the verifier
        isfound = true
    end

    if isfound
        println("Valid CLF: terminated")
        return BARRIER_FOUND, gen_prob
    end

    if !issuccess
        println("No valid CLF: terminated")
        return BARRIER_INFEASIBLE, gen_prob
    end

    if iter ≥ iter_max
        println("\nMax iter exceeded: ", iter)
        return MAX_ITER_REACHED, gen_prob
    end

    error("Something bad")
end