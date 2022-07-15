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
    neg::PointSet{N,M}
    unsafe::PointSet{N,M}
    pos::PointSet{N,M}
end

Witness{N,M}() where {N,M} = Witness(ntuple(k -> PointSet{N,M}(), Val(3))...)

# Recorder

snapshot(::Nothing, ::Any) = nothing

struct TraceRecorder{N,M}
    wit_list::Vector{Witness{N,M}}
end

TraceRecorder{N,M}() where {N,M} = TraceRecorder(Witness{N,M}[])
TraceRecorder(::Learner{N,M}) where {N,M} = TraceRecorder{N,M}()

function snapshot(rec::TraceRecorder, wit::Witness)
    neg_ = PointSet(copy.(wit.neg.points_list))
    unsafe_ = PointSet(copy.(wit.unsafe.points_list))
    pos_ = PointSet(copy.(wit.pos.points_list))
    push!(rec.wit_list, Witness(neg_, unsafe_, pos_))
end

## Learn Barrier
function learn_lyapunov!(
        lear::Learner{N,M}, iter_max, solver_sep, solver_verif;
        do_print=true, rec=nothing
    ) where {N,M}
    ϵ = lear.ϵ
    βmax = lear.params[:βmax]
    xmax = lear.params[:xmax]
    tol_dom = lear.params[:tol_dom]

    wit = Witness{N,M}()
    temp_pos_set = PointSet{N,M}()

    for (loc, points) in enumerate(lear.iset.points_list)
        for point in points
            add_point!(wit.neg, loc, point)
        end
    end

    mpf = MultiPolyFunc{N,M}()
    iter = 0
    loc_stack = collect(1:M)
    
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
        neg_points = wit.neg.points_list[loc]

        # Sep unsafe
        for point in wit.unsafe.points_list[loc]
            af, r = compute_af(neg_points, point, ϵ, βmax, solver_sep)
            if r < 0
                println(string("Satisfiability radius too small: ", r))
                return RADIUS_TOO_SMALL, mpf, wit
            end
            add_af!(mpf, loc, af)
        end

        empty!(temp_pos_set, loc)

        while !isempty(wit.pos.points_list[loc])
            point = pop_point!(wit.pos, loc)
            af, r = compute_af(neg_points, point, ϵ, βmax, solver_sep)
            if r < 0
                println(string("|--- radius: ", r))
                for piece in lear.sys.pieces
                    loc != piece.loc1 && continue
                    !_neg(piece.pf_dom, point, tol_dom) && continue
                    loc2 = piece.loc2
                    add_point!(wit.neg, loc2, piece.A*point + piece.b)
                    loc2 ∉ loc_stack && push!(loc_stack, loc2)
                end
                break
            end
            add_af!(mpf, loc, af)
            add_point!(temp_pos_set, loc, point)
        end

        for point in temp_pos_set.points_list[loc]
            add_point!(wit.pos, loc, point)
        end

        !isempty(loc_stack) && continue

        # Verifier
        verif = Verifier(lear.mpf_safe, lear.mpf_inv, mpf, lear.sys, xmax)
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(verif, lear.δ, solver_verif)
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            add_point!(wit.unsafe, loc, x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(verif, lear.δ, solver_verif)
        if obj > 0
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            add_point!(wit.pos, loc, x)
            push!(loc_stack, loc)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
    end
    println("Valid CLF: terminated")
    return BARRIER_FOUND, mpf, wit
end