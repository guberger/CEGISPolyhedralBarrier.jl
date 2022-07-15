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

# Witness
struct Witness{N,M}
    inside::PointSet{N,M}
    image::PointSet{N,M}
    outside::PointSet{N,M}
    unknown::PointSet{N,M}
end

Witness{N,M}() where {N,M} = Witness(ntuple(k -> PointSet{N,M}(), Val(4))...)

# Recorder

snapshot(::Nothing, ::Any) = nothing

struct TraceRecorder{N,M}
    wit_list::Vector{Witness{N,M}}
end

TraceRecorder{N,M}() where {N,M} = TraceRecorder(Witness{N,M}[])
TraceRecorder(::Learner{N,M}) where {N,M} = TraceRecorder{N,M}()

function snapshot(rec::TraceRecorder, wit::Witness)
    inside_ = PointSet(copy.(wit.inside.points_list))
    image_ = PointSet(copy.(wit.image.points_list))
    outside_ = PointSet(copy.(wit.outside.points_list))
    unknown_ = PointSet(copy.(wit.unknown.points_list))
    push!(rec.wit_list, Witness(inside_, image_, outside_, unknown_))
end

# Add new point to invariant set
function _add_safe_point!(wit, loc_stack, sys, loc, point, tol_dom)
    add_point!(wit.inside, loc, point)
    loc ∉ loc_stack && push!(loc_stack, loc)
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !_prox(piece.pf_dom, point, tol_dom) && continue
        loc2 = piece.loc2
        add_point!(wit.image, loc2, piece.A*point + piece.b)
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
        lear::Learner{N,M}, iter_max, solver_sep, solver_verif;
        do_print=true, rec=nothing
    ) where {N,M}
    sys = lear.sys
    ϵ = lear.ϵ
    δ = lear.δ
    βmax = lear.params[:βmax]
    xmax = lear.params[:xmax]
    tol_dom = lear.params[:tol_dom]

    wit = Witness{N,M}()
    temp_unknown_set = PointSet{N,M}()
    loc_stack = Int[]

    for (loc, points) in enumerate(lear.iset.points_list)
        for point in points
            _add_safe_point!(wit, loc_stack, sys, loc, point, tol_dom)
        end
    end

    mpf = MultiPolyFunc{N,M}()
    iter = 0
    
    while !isempty(loc_stack)        
        iter += 1
        do_print && println("Iter: ", iter)
        if iter > iter_max
            println(string("Max iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, wit
        end

        snapshot(rec, wit)

        loc = pop!(loc_stack)
        empty!(mpf, loc)
        inside_points = wit.inside.points_list[loc]
        image_points = wit.image.points_list[loc]

        # Sep outside
        for point in wit.outside.points_list[loc]
            af, r = compute_af(image_points, point, ϵ, βmax, solver_sep)
            if r < 0
                println(string("Satisfiability radius too small: ", r))
                return RADIUS_TOO_SMALL, mpf, wit
            end
            add_af!(mpf, loc, af)
            af, r = compute_af(inside_points, point, 0, βmax, solver_sep)
            if r < 0
                println(string("Satisfiability radius too small: ", r))
                return RADIUS_TOO_SMALL, mpf, wit
            end
        end

        empty!(temp_unknown_set, loc)

        while !isempty(wit.unknown.points_list[loc])
            point = pop_point!(wit.unknown, loc)
            af, r = compute_af(image_points, point, ϵ, βmax, solver_sep)
            if r < 0
                println(string("|--- radius: ", r))
                _add_safe_point!(wit, loc_stack, sys, loc, point, tol_dom)
                break
            end
            add_af!(mpf, loc, af)
            af, r = compute_af(inside_points, point, 0, βmax, solver_sep)
            if r < 0
                println(string("|--- radius: ", r))
                _add_safe_point!(wit, loc_stack, sys, loc, point, tol_dom)
                break
            end
            add_point!(temp_unknown_set, loc, point)
        end

        for point in temp_unknown_set.points_list[loc]
            add_point!(wit.unknown, loc, point)
        end

        !isempty(loc_stack) && continue

        # Verifier
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(
            sys, lear.mpf_safe, lear.mpf_inv, mpf, xmax, -δ, solver_verif
        )
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            add_point!(wit.outside, loc, x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(
            sys, lear.mpf_safe, lear.mpf_inv, mpf, xmax, -δ, solver_verif
        )
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            if _is_outside(sys, lear.mpf_safe, loc, x, ϵ, tol_dom)
                @assert false
                add_point!(wit.outside, loc, x)
            else
                add_point!(wit.unknown, loc, x)
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