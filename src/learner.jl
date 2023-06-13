## Learner

@enum StatusCode begin
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    MAX_ITER_REACHED = 3
end

function _reset_crossing_problem!(prob)
    empty!(prob.afs_inside)
    empty!(prob.afs_inside_margin)
    empty!(prob.afs_outside)
    empty!(prob.afs_outside_margin)
end

struct VerifierSafeProblem
    N::Int
    sys::System
    mpf_safe::MultiPolyFunc
    mpf_inv::MultiPolyFunc
    mpf_BF::MultiPolyFunc
end

function find_counterexample(prob::VerifierSafeProblem, xmax, solver)
    cross_prob = empty_cross_problem(prob.N)
    xopt::Vector{Float64} = fill(NaN, prob.N)
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in prob.sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        _reset_crossing_problem!(cross_prob)
        append!(cross_prob.afs_inside, pf_dom.afs)
        append!(cross_prob.afs_inside, prob.mpf_inv.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_safe.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_BF.pfs[piece.loc1].afs)
        for af in prob.mpf_safe.pfs[piece.loc2].afs
            empty!(cross_prob.afs_outside)
            push!(cross_prob.afs_outside, af)
            x, r, flag = find_crosser(cross_prob, A, b, xmax, solver
            )
            if flag && r > ropt
                xopt = x
                ropt = r
                locopt = piece.loc1
            end
        end
    end
    return xopt, ropt, locopt
end

struct VerifierBFProblem
    N::Int
    sys::System
    mpf_safe::MultiPolyFunc
    mpf_inv::MultiPolyFunc
    mpf_BF::MultiPolyFunc
end

function find_counterexample(prob::VerifierBFProblem, xmax, solver)
    cross_prob = empty_cross_problem(prob.N)
    xopt::Vector{Float64} = fill(NaN, prob.N)
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in prob.sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        _reset_crossing_problem!(cross_prob)
        append!(cross_prob.afs_inside, pf_dom.afs)
        append!(cross_prob.afs_inside, prob.mpf_inv.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_safe.pfs[piece.loc1].afs)
        append!(cross_prob.afs_inside_margin, prob.mpf_BF.pfs[piece.loc1].afs)
        for af in prob.mpf_BF.pfs[piece.loc2].afs
            empty!(cross_prob.afs_outside_margin)
            push!(cross_prob.afs_outside_margin, af)
            x, r, flag = find_crosser(cross_prob, A, b, xmax, solver
            )
            if flag && r > ropt
                xopt = x
                ropt = r
                locopt = piece.loc1
            end
        end
    end
    return xopt, ropt, locopt
end

struct BarrierProblem
    N::Int
    M::Int
    sys::System,
    mpf_safe::MultiPolyFunc
    mpf_inv::MultiPolyFunc
    mgrid_init::Vector{Vector{Vector{Float64}}},
    ϵ::Float64
end

# :tol_dom=1e-8,
# :βmax=1e3,
# :xmax=1e3

struct LearnerState
    sep_prob::SeparationProblem
    verif_safe_prob::VerifierSafeProblem
    verif_bf_prob::VerifierBFProblem
    mgrid_inside::Vector{Vector{Vector{Float64}}}
    mgrid_image::Vector{Vector{Vector{Float64}}}
    mgraph_unknown::Vector{Vector{Vector{Float64}}}
    mgraph_unknown_new::Vector{Vector{Vector{Float64}}}
    mgrid_outside::Vector{Vector{Vector{Float64}}}
    mgrid_outside_new::Vector{Vector{Vector{Float64}}}
    loc_stack::Vector{Int}
end
state_init(M) = Cluster(ntuple(i -> empty_multigrid(M), Val(6))..., Int[])

# Add new point to invariant set
function add_loc!(loc_stack, loc)
    if loc ∉ state.loc_stack
        push!(state.loc_stack, loc)
    end
end

function add_safe_point!(state::State, sys, loc, point, tol_dom)
    push!(state.mgrid_inside[loc], point)
    add_loc!(state.loc_stack, loc)
    for piece in sys.pieces
        if loc == piece.loc1 && isless_approx(piece.pf_dom, point, 0, tol_dom)
            push!(state.mgrid_image[piece.loc2], piece.A*point + piece.b)
            add_loc!(state.loc_stack, piece.loc2)
        end
    end
end

function _is_outside(sys, mpf_safe, loc, point, ϵ, tol_dom)
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !_prox(piece.pf_dom, point, tol_dom) && continue
        loc2 = piece.loc2
        point2 = piece.A*point + piece.b
        if !_prox(mpf_safe.pfs[loc2], point2, -2*ϵ)
            return true
        end
    end
    return false
end

function update_mpf!(mpf, prob, state, ϵ, βmax, solver)
    while !isempty(state.loc_stack)
        for point in state.mgrid_outside_new[loc]
            af, r = compute_af(
                grid_inside, grid_image, point, βmax, N, solver_sep
            )
            if r < ϵ
                return false
            end
            push!(mpf.pfs[loc].afs, af)
        end
end

## Learn Barrier
function find_barrier(prob::BarrierProblem, iter_max, solver_sep, solver_verif;
                      tol_dom=1e-8, βmax=1e3, xmax=1e3,
                      do_print=true, callback_fcn=(args...) -> nothing)
    sys, mpf_safe, mpf_inv = prob.sys, prob.mpf_safe, prob.mpf_inv
    state = state_init(prob.M)
    for loc = 1:prob.M
        for point in mgrid_init[loc]
            add_safe_point!(state, sys, loc, point, tol_dom)
        end
    end
    mpf = empty_mpf(prob.M)
    prob = SeparationProblem(prob.N, )
    iter = 0
    
    while !isempty(state.loc_stack)        
        iter += 1
        do_print && print("Iter: ", iter)
        if iter > iter_max
            println(string("\nMax iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, state
        end

        callback_fcn(iter, mpf, state)

        loc = pop!(state.loc_stack)
        do_print && println(" - loc:", loc)
        empty!(mpf, loc)
        grid_inside = state.mgrid_inside[loc]
        grid_image = state.mgrid_image[loc]


        # Sep outside
        for point in state.mgrid_outside[loc]
            af, r = compute_af(
                grid_inside, grid_image, point, βmax, N, solver_sep
            )
            if r < ϵ
                println(string("Radius too small: ", r))
                return BARRIER_INFEASIBLE, mpf, state
            end
            push!(mpf.pfs[loc].afs, af)
        end

        empty!(list_unknown_temp)

        while !isempty(state.mgraph_unknown[loc])
            point = pop!(state.mgraph_unknown[loc])
            af, r = compute_af(
                grid_inside, grid_image, point, βmax, N, solver_sep
            )
            if r < ϵ
                do_print && println(string("|--- radius: ", r))
                _add_safe_point!(
                    state.mgrid_inside, state.mgrid_image, loc_stack,
                    _sys, loc, point, tol_dom
                )
                # uncomment two lines below for approach of the paper
                # empty!(state.mgraph_unknown[loc])
                # empty!(list_unknown_temp)
                break
            end
            push!(mpf.pfs[loc].afs, af)
            push!(list_unknown_temp, point)
        end

        for point in list_unknown_temp
            push!(state.mgraph_unknown[loc], point)
        end

        !isempty(loc_stack) && continue

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

    println("Valid CLF: terminated")
    return BARRIER_FOUND, mpf, state
end