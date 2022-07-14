## Learner

@enum StatusCode begin
    NOT_SOLVED = 0
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    RADIUS_TOO_SMALL = 3
    MAX_ITER_REACHED = 4
end

## Learner

struct Learner{N,M}
    sys::System{N}
    mpf_safe::MultiPolyFunc{N,M}
    mpf_inv::MultiPolyFunc{N,M}
    iset::PointSet{N,M}
    ϵ::Float64
    δ::Float64
    params::Dict{Symbol,Float64}
end

function Learner(sys, mpf_safe, mpf_inv, isetϵ, ϵ, δ)
    params = Dict([
        :tol_dom => 1e-8,
        :βmax => 1e3,
        :xmax => 1e3
    ])
    return Learner(sys, mpf_safe, mpf_inv, isetϵ, ϵ, δ, params)
end

_setsafe!(D, k, v) = (@assert haskey(D, k); D[k] = v)
set_param!(lear::Learner, s::Symbol, v) = _setsafe!(lear.params, s, v)

## Learn Barrier
function learn_lyapunov!(
        lear::Learner{N,M}, iter_max, solver_sep, solver_verif; do_print=true
    ) where {N,M}
    βmax = lear.params[:βmax]
    xmax = lear.params[:xmax]
    tol_dom = lear.params[:tol_dom]

    seps = ntuple(loc -> Separator{N}(lear.ϵ, βmax), Val(M))
    for (sep, points) in zip(seps, lear.iset.points_list)
        for point in points
            add_soft_evid!(sep, point)
        end
    end

    unsafe_set = PointSet{N,M}()
    pos_set = PointSet{N,M}()
    temp_pos_set = PointSet{N,M}()

    mpf = MultiPolyFunc{N,M}()
    iter = 0
    loc_stack = collect(1:M)
    
    while !isempty(loc_stack)        
        iter += 1
        do_print && println("Iter: ", iter)
        if iter > iter_max
            println(string("Max iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, seps, iter
        end

        loc = pop!(loc_stack)
        sep = seps[loc]
        empty!(mpf, loc)

        # Sep unsafe
        for point in unsafe_set.points_list[loc]
            af, r = compute_af(sep, point, solver_sep)
            if r < 0
                println(string("Satisfiability radius too small: ", r))
                return RADIUS_TOO_SMALL, mpf, seps, iter
            end
            add_af!(mpf, loc, af)
        end

        empty!(temp_pos_set, loc)

        while !isempty(pos_set.points_list[loc])
            point = pop_point!(pos_set, loc)
            af, r = compute_af(sep, point, solver_sep)
            if r < 0
                println(string("|--- radius: ", r))
                add_soft_evid!(sep, point)
                for piece in lear.sys.pieces
                    loc != piece.loc1 && continue
                    !_neg(piece.pf_dom, point, tol_dom) && continue
                    loc2 = piece.loc2
                    add_hard_evid!(seps[loc2], piece.A*point + piece.b)
                    loc2 ∉ loc_stack && push!(loc_stack, loc2)
                end
                break
            end
            add_af!(mpf, loc, af)
            add_point!(temp_pos_set, loc, point)
        end

        for point in temp_pos_set.points_list[loc]
            add_point!(pos_set, loc, point)
        end

        !isempty(loc_stack) && continue
        @assert isempty(loc_stack)

        # Verifier
        verif = Verifier(lear.mpf_safe, lear.mpf_inv, mpf, lear.sys, xmax)
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(verif, lear.δ, solver_verif)
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            add_point!(unsafe_set, loc, x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(verif, lear.δ, solver_verif)
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            add_point!(pos_set, loc, x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
    end
    println("Valid CLF: terminated")
    return BARRIER_FOUND, mpf, seps, iter
end