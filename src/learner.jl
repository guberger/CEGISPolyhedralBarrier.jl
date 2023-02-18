## Learner

@enum StatusCode begin
    NOT_SOLVED = 0
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    MAX_ITER_REACHED = 3
end

# :tol_dom=1e-8,
# :βmax=1e3,
# :xmax=1e3

struct Witness
    mlist_inside::Vector{Vector{Vector{Float64}}}
    mlist_image::Vector{Vector{Vector{Float64}}}
    mlist_unknown::Vector{Vector{Vector{Float64}}}
    mlist_outside::Vector{Vector{Vector{Float64}}}
end

# Add new point to invariant set
function _add_safe_point!(
        mlist_inside, mlist_image, loc_stack, sys, loc, point, tol_dom
    )
    push!(mlist_inside[loc], point)
    loc ∉ loc_stack && push!(loc_stack, loc)
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !_prox(piece.pf_dom, point, tol_dom) && continue
        loc2 = piece.loc2
        push!(mlist_image[loc2], piece.A*point + piece.b)
        loc2 ∉ loc_stack && push!(loc_stack, loc2)
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

## Learn Barrier
function learn_lyapunov!(
        sys::System,
        mpf_safe::MultiPolyFunc, mpf_inv::MultiPolyFunc,
        mlist_init::AbstractVector{<:AbstractVector{<:AbstractVector}},
        ϵ, δ, iter_max, M, N, solver_sep, solver_verif;
        tol_dom=1e-8, βmax=1e3, xmax=1e3,
        do_print=true, callback_fcn=(args...) -> nothing
    )
    _mpf_safe = MultiPolyFunc([
        PolyFunc([
            AffForm(Vector{Float64}(af.a), Float64(af.β))
            for af in mpf_safe.pfs[loc].afs
        ]) for loc = 1:M
    ])
    _mpf_inv = MultiPolyFunc([
        PolyFunc([
            AffForm(Vector{Float64}(af.a), Float64(af.β))
            for af in mpf_inv.pfs[loc].afs
        ]) for loc = 1:M
    ])

    wit = Witness(ntuple(i -> [Vector{Float64}[] for loc = 1:M], Val(4))...)
    loc_stack = Int[]
    for loc = 1:M
        for point in mlist_init[loc]
            _add_safe_point!(
                wit.mlist_inside, wit.mlist_image, loc_stack,
                sys, loc, point, tol_dom
            )
        end
    end

    mpf = MultiPolyFunc([
        PolyFunc(AffForm{Vector{Float64},Float64}[]) for loc = 1:M
    ])
    iter = 0
    list_unknown_temp = Vector{Float64}[]
    
    while !isempty(loc_stack)        
        iter += 1
        do_print && print("Iter: ", iter)
        if iter > iter_max
            println(string("\nMax iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, wit
        end

        callback_fcn(iter, mpf, wit)

        loc = pop!(loc_stack)
        do_print && println(" - loc:", loc)
        empty!(mpf, loc)
        points_inside = wit.mlist_inside[loc]
        points_image = wit.mlist_image[loc]

        # Sep outside
        for point in wit.mlist_outside[loc]
            af, r = compute_af(
                points_inside, points_image, point, ϵ, βmax, N, solver_sep
            )
            if r < 0
                println(string("Radius too small: ", r))
                return BARRIER_INFEASIBLE, mpf, wit
            end
            push!(mpf.pfs[loc].afs, af)
        end

        empty!(list_unknown_temp)

        while !isempty(wit.mlist_unknown[loc])
            point = pop!(wit.mlist_unknown[loc])
            af, r = compute_af(
                points_inside, points_image, point, ϵ, βmax, N, solver_sep
            )
            if r < 0
                println(string("|--- radius: ", r))
                _add_safe_point!(
                    wit.mlist_inside, wit.mlist_image, loc_stack,
                    sys, loc, point, tol_dom
                )
                break
            end
            push!(mpf.pfs[loc].afs, af)
            push!(list_unknown_temp, point)
        end

        for point in list_unknown_temp
            push!(wit.mlist_unknown[loc], point)
        end

        !isempty(loc_stack) && continue

        # Verifier
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(
            sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
        )
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            push!(wit.mlist_outside[loc], x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(
            sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
        )
        if obj > 0
            do_print && print("CE found: ", x, ", ", loc, ", ", obj)
            if _is_outside(sys, _mpf_safe, loc, x, ϵ, tol_dom)
                println(" (outside)")
                push!(wit.mlist_outside[loc], x)
            else
                println(" (unknown)")
                push!(wit.mlist_unknown[loc], x)
            end
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
    end

    println("Valid CLF: terminated")
    return BARRIER_FOUND, mpf, wit
end