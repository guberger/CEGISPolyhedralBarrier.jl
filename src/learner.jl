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

const _F_AffForm = AffForm{Vector{Float64},Float64}
const _F_PolyFunc = PolyFunc{_F_AffForm}
const _F_MultiPolyFunc = MultiPolyFunc{_F_PolyFunc}
const _F_Piece = Piece{_F_PolyFunc,Matrix{Float64},Vector{Float64}}
const _F_System = System{_F_Piece}

function _float_af(af::AffForm)::_F_AffForm
    return AffForm(Vector{Float64}(af.a), Float64(af.β))
end

function _float_pf(pf::PolyFunc)::_F_PolyFunc
    return PolyFunc([_float_af(af) for af in pf.afs])
end

function _float_mpf(mpf::MultiPolyFunc)::_F_MultiPolyFunc
    return MultiPolyFunc([_float_pf(pf) for pf in mpf.pfs])
end

function _float_piece(piece::Piece)::_F_Piece
    return Piece(
        _float_pf(piece.pf_dom), piece.loc1,
        Matrix{Float64}(piece.A), Vector{Float64}(piece.b), piece.loc2
    )
end

function _float_sys(sys::System)::_F_System
    return System([_float_piece(piece) for piece in sys.pieces])
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
    _mpf_safe = _float_mpf(mpf_safe)
    _mpf_inv = _float_mpf(mpf_inv)
    _sys = _float_sys(sys)

    wit = Witness(ntuple(i -> [Vector{Float64}[] for loc = 1:M], Val(4))...)
    loc_stack = Int[]
    for loc = 1:M
        for point in mlist_init[loc]
            _add_safe_point!(
                wit.mlist_inside, wit.mlist_image, loc_stack,
                _sys, loc, point, tol_dom
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
                points_inside, points_image, point, βmax, N, solver_sep
            )
            if r < ϵ
                println(string("Radius too small: ", r))
                return BARRIER_INFEASIBLE, mpf, wit
            end
            push!(mpf.pfs[loc].afs, af)
        end

        empty!(list_unknown_temp)

        while !isempty(wit.mlist_unknown[loc])
            point = pop!(wit.mlist_unknown[loc])
            af, r = compute_af(
                points_inside, points_image, point, βmax, N, solver_sep
            )
            if r < ϵ
                do_print && println(string("|--- radius: ", r))
                _add_safe_point!(
                    wit.mlist_inside, wit.mlist_image, loc_stack,
                    _sys, loc, point, tol_dom
                )
                # uncomment two lines below for approach of the paper
                # empty!(wit.mlist_unknown[loc])
                # empty!(list_unknown_temp)
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
            _sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
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
            _sys, _mpf_safe, _mpf_inv, mpf, xmax, -δ, N, solver_verif
        )
        if obj > 0
            do_print && print("CE found: ", x, ", ", loc, ", ", obj)
            if _is_outside(_sys, _mpf_safe, loc, x, ϵ, tol_dom)
                do_print && println(" (outside)")
                push!(wit.mlist_outside[loc], x)
            else
                do_print && println(" (unknown)")
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