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
end

# :tol=1e-8,
# :βmax=1e3,
# :xmax=1e3

function build_generator(prob, βmax, isint)
    gfs_safe = copy(prob.gfs_safe)
    for (i, gf) in enumerate(gfs_safe)
        na = norm(gf.af.a, Inf)
        @assert na > 1e-5
        gfs_safe[i] = GenForm(gf.loc, AffForm(gf.af.a / na, gf.af.β / na))
    end
    return GeneratorProblem(
        prob.N,
        GenForm[], # gfs
        gfs_safe, # gfs_safe
        Int[], # indices_new
        copy(prob.states_init), # states_inside
        State[], # states_image
        Edge[], # edges_unknown
        Edge[], # edges_unknown_new
        [Vector{Float64}[] for i = 1:3]..., # all xs_...
        prob.ϵ, βmax, isint
    )
end

function build_verifier(prob, xmax, isint)
    gfs_inv = copy(prob.gfs_inv)
    for (i, gf) in enumerate(gfs_inv)
        na = norm(gf.af.a, Inf)
        @assert na > 1e-5
        gfs_inv[i] = GenForm(gf.loc, AffForm(gf.af.a / na, gf.af.β / na))
    end
    return VerifierProblem(
        prob.N,
        prob.pieces,
        gfs_inv,
        GenForm[], # gfs_new
        GenForm[], # gfs_inside
        GenForm[], # gfs_outside
        Int[], # indices_inside_new
        Int[], # indices_outside_new
        Witness[], # witnesses
        Tuple{Int,Int,Int}[], # keys_todo
        [AffForm[] for i = 1:2]..., # all afs_...
        xmax, isint
    )
end

## Learn Barrier
function find_barrier(prob::BarrierProblem,
                      iter_max, solver; # LP solver
                      δ=1e-8, βmax=1e3, xmax=1e3,
                      isint_gen=false, isint_verif=false)
    gen_prob = build_generator(prob, βmax, isint_gen)
    verif_prob = build_verifier(prob, xmax, isint_verif)
    iter = 0
    isfound = false
    issuccess = true
    
    while iter < iter_max && !isfound && issuccess
        iter += 1
        @printf("iter %4d: ", iter)
        time_start = time()

        # Generation part
        @printf("xs: [%d %d %d %d] ",
                length(gen_prob.states_inside),
                length(gen_prob.states_image),
                length(gen_prob.edges_unknown),
                length(gen_prob.edges_unknown_new))
        isreset, issuccess = update_generator!(gen_prob, solver)
        issuccess || break
        @assert isempty(gen_prob.edges_unknown_new)
        @printf("reset: %d ", isreset)

        # Verification part
        @printf("gfs: %d %d ",
                length(gen_prob.gfs),
                length(gen_prob.indices_new))
        empty!(verif_prob.gfs_new)
        for i in gen_prob.indices_new
            push!(verif_prob.gfs_new, gen_prob.gfs[i])
        end
        kopt, ropt = update_verifier!(verif_prob, isreset, solver)
        @printf("ropt: %f ", ropt)
        if ropt > -δ
            witness = verif_prob.witnesses[kopt]
            state = witness.state
            piece = verif_prob.pieces[witness.q]
            state_post = State(piece.loc_dst, piece.A*state.x + piece.b)
            edge = Edge(state, state_post)
            push!(gen_prob.edges_unknown_new, edge)
        else
            isfound = true
        end

        time_elapsed = time() - time_start
        @printf("time: %f\n", time_elapsed)
    end

    if isfound
        @printf("Valid BF: terminated\n")
        return BARRIER_FOUND, gen_prob
    end
    if !issuccess
        @printf("No valid BF: terminated\n")
        return BARRIER_INFEASIBLE, gen_prob
    end
    if iter ≥ iter_max
        @printf("\nMax iter exceeded: %d\n", iter)
        return MAX_ITER_REACHED, gen_prob
    end
    error("Something bad")
end