## Learner

@enum StatusCode begin
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    MAX_ITER_REACHED = 3
end

struct BarrierProblem
    N::Int
    pieces::Vector{Piece}
    gfs_safe::Vector{GenForm}
    gfs_inv::Vector{GenForm}
    states_init::Vector{State}
    ϵ::Float64
end

# :tol_dom=1e-8,
# :βmax=1e3,
# :xmax=1e3

function build_generator(prob)
    return GeneratorProblem(
        prob.N,
        GenForm[], # gfs
        int[], # indices_new
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

function build_verifier(prob)
    return VerifierProblem(
        prob.N,
        prob.pieces,
        GenForm[], # gfs_bf
        prob.gfs_safe,
        prob.gfs_inv,
        Dict{CexKey,CexVal}(), # cexs
        CexKey[], # keys_todo
        [AffForm[] for i = 1:4]..., # all afs_...
    )
end

function reset_verifier_safe!(prob::VerifierProblem)
    empty!(prob.cexs)
    empty!(prob.keys_todo)
    add_gfs_keys!(prob, prob.gfs_safe, eachindex(prob.gfs_safe))
end

function reset_verifier_cont!(prob::VerifierProblem)
    empty!(prob.cexs)
    empty!(prob.keys_todo)
    add_gfs_keys!(prob, prob.gfs_bf, eachindex(prob.gfs_bf))
end

function copy_filter!(gfs_new, gfs, indices)
    
end

## Learn Barrier
function find_barrier(prob::BarrierProblem,
                      iter_max, solver; # LP solver
                      βmax=1e3, xmax=1e3,
                      do_print=true, callback_fcn=(args...) -> nothing)
    
    gen_prob = build_generator(prob)
    verif_safe_prob = build_verifier(prob)
    reset_verifier_safe!(verif_safe_prob)
    verif_cont_prob = build_verifier(prob)
    reset_verifier_cont!(verif_cont_prob)
    gfs_new = GenForm[]

    iter = 0
    isfound = false
    issuccess = true
    
    while iter < iter_max && !isfound && issuccess
        iter += 1
        do_print && print("Iter: ", iter)

        callback_fcn(iter, gen_prob)

        # Generation part
        isreset, issuccess = update_generator!(prob, βmax, solver)

        if !issuccess
            break
        end

        # Verification safe part
        if isreset
            reset_verifier_safe!(verif_safe_prob.keys_todo)
        end
        update_cexs_safe!(verif_safe_prob, xmax, solver)


        # Verifier
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(
            _sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
        )
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            push!(state.mgrid_outside[loc], x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(
            _sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
        )
        if obj > 0
            do_print && print("CE found: ", x, ", ", loc, ", ", obj)
            if _is_outside(_sys, _mpf_safe, loc, x, ϵ, tol_dom)
                do_print && println(" (outside)")
                push!(state.mgrid_outside[loc], x)
            else
                do_print && println(" (unknown)")
                push!(state.mgraph_unknown[loc], x)
            end
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
    end

    if iter > iter_max
        println(string("\nMax iter exceeded: ", iter))
        return MAX_ITER_REACHED, mpf, state
    end

    println("Valid CLF: terminated")
    return BARRIER_FOUND, mpf, state
end