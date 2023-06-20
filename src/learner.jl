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
        copy(prob.gfs_safe), # gfs_safe
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

function build_verifier(prob)
    return VerifierProblem(
        prob.N,
        prob.pieces,
        copy(prob.gfs_inv),
        GenForm[], # gfs_bf
        GenForm[], # gfs_out
        Dict{CexKey,CexVal}(), # cexs
        CexKey[], # keys_todo
        [AffForm[] for i = 1:2]... # all afs_...
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

struct PrintRecorder
    ninside::Vector{Int}
    nimage::Vector{Int}
    nunknown::Vector{Int}
    noutside::Vector{Int}
    ngf::Vector{Int}
    ngf_new::Vector{Int}
    gen_res::Vector{Char}
    nkey::Vector{Int}
    nkey_todo::Vector{Int}
    rs::Vector{Float64}
end

function init_recorder()
    return PrintRecorder(
        [Int[] for i = 1:6]..., Char[], [Int[] for i = 1:2]..., Float64[]
    )
end

function reset_recorder!(rec::PrintRecorder)
    empty!(rec.ninside); empty!(rec.nimage);
    empty!(rec.nunknown); empty!(rec.noutside);
    empty!(rec.ngf); empty!(rec.ngf_new);
    empty!(rec.gen_res);
    empty!(rec.nkey); empty!(rec.nkey_todo); empty!(rec.rs)
end

function update_recorder!(rec::PrintRecorder, prob::GeneratorProblem)
    push!(rec.ninside, length(prob.states_inside))
    push!(rec.nimage, length(prob.states_image))
    push!(rec.nunknown, length(prob.links_unknown))
    push!(rec.noutside, length(prob.states_outside))
    push!(rec.ngf, length(prob.gfs))
    push!(rec.ngf_new, length(prob.indices_new))
end

function update_recorder!(rec::PrintRecorder, prob::VerifierProblem)
    push!(rec.nkey, length(prob.cexs))
    push!(rec.nkey_todo, length(prob.keys_todo))
end

function print_algo_state(iter::Int, rec::PrintRecorder)
    print(
        "Iter: ", iter, "\n",
        " - ninside: ", rec.ninside, "\n",
        " - nimage: ", rec.nimage, "\n",
        " - nunknown: ", rec.nunknown, "\n",
        " - noutside: ", rec.noutside, "\n",
        " - ngf: ", rec.ngf, "\n",
        " - ngf_new: ", rec.ngf_new, "\n",
        " - gen_res: ", rec.gen_res, "\n",
        " - nkey: ", rec.nkey, "\n",
        " - nkey_todo: ", rec.nkey_todo, "\n"
    )
    if !isempty(rec.rs)
        print(" - rs: ", minimum(rec.rs), " -- ", maximum(rec.rs), "\n")
    end
end

## Learn Barrier
function find_barrier(prob::BarrierProblem,
                      iter_max, solver; # LP solver
                      βmax=1e3, xmax=1e3,
                      print_period::Int=1)
    
    gen_prob = build_generator(prob)
    verif_prob = build_verifier(prob)
    rec = init_recorder()

    iter = 0
    isfound = false
    issuccess = true
    
    while iter < iter_max && !isfound && issuccess
        iter += 1

        # Generation part
        isreset, issuccess = update_generator!(gen_prob, βmax, solver)
        @assert isempty(gen_prob.links_unknown_new)
        @assert isempty(gen_prob.states_outside_new)

        update_recorder!(rec, gen_prob)

        if !issuccess
            break
        end

        # Verification part
        copy!(verif_prob.gfs_bf, gen_prob.gfs)
        copy!(verif_prob.gfs_out, gen_prob.gfs)
        if isreset
            reset_verifier!(verif_prob)
            push!(rec.gen_res, 'R')
        else
            empty!(verif_prob.keys_todo)
            add_keys_bf_infeasible!(verif_prob, gen_prob.indices_new)
            add_keys_out_new!(verif_prob, gen_prob.indices_new)
            push!(rec.gen_res, 'E')
        end

        update_cexs!(verif_prob, xmax, solver)
        update_recorder!(rec, verif_prob)

        key, r = find_cex_max(verif_prob.cexs)
        push!(rec.rs, r)
        
        if r > -prob.δ
            state = verif_prob.cexs[key].state
            piece = verif_prob.pieces[key.q]
            state_post = State(piece.loc_dst, piece.A*state.x + piece.b)
            link = Link(state, state_post)
            push!(gen_prob.links_unknown_new, link)
        else
            isfound = true
        end

        if print_period > 0 && mod(iter, print_period) == 0
            print_algo_state(iter, rec)
            reset_recorder!(rec)
        end
    end

    print_algo_state(iter, rec)

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