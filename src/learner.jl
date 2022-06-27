## Learner

@enum StatusCode begin
    NOT_SOLVED = 0
    BARRIER_FOUND = 1
    BARRIER_INFEASIBLE = 2
    RADIUS_TOO_SMALL = 3
    MAX_ITER_REACHED = 4
end

## Learner

struct Learner
    nvar::Int
    nloc::Int
    sys::System
    iset::InitSet
    uset::UnsafeSet
    ϵ::Float64
    δ::Float64
    tols::Dict{Symbol,Float64}
    params::Dict{Symbol,Float64}
end

function Learner(
        nvar::Int, nloc::Int,
        sys::System, iset::InitSet, uset::UnsafeSet,
        ϵ::Float64, δ::Float64
    )
    tols = Dict([
        :rad => eps(1.0),
        :pos => -eps(1.0),
        :lie => -eps(1.0)
    ])
    params = Dict([
        :offmax => 1e3,
        :rmax_gen => 1e4,
        :bigM => 1e3,
        :xmax => 1e3,
        :rmax_verif => 1e2
    ])
    return Learner(nvar, nloc, sys, iset, uset, ϵ, δ, tols, params)
end

set_tol!(lear::Learner, s::Symbol, v::Float64) = (lear.tols[s] = v)
set_param!(lear::Learner, s::Symbol, v::Float64) = (lear.params[s] = v)

## Learn Barrier

function _add_evidences_pos!(gen, state)
    i = add_af!(gen, state.loc)
    add_evidence!(gen, PosEvidence(state.loc, i, state.point))
end

function _add_evidences_neg!(gen, state)
    add_evidence!(gen, NegEvidence(state.loc, state.point))
end

function _add_evidences_lie!(gen, sys, state)
    point1 = state.point
    loc1 = state.loc
    i1 = add_af!(gen, loc1)
    for piece in sys.pieces
        !(loc1 == piece.loc1 && point1 ∈ piece.domain) && continue
        point2 = piece.A*point1 + piece.b
        loc2 = piece.loc2
        nA = opnorm(piece.A, Inf)
        add_evidence!(gen, LieEvidence(loc1, i1, point1, loc2, point2, nA))
    end
end

function _add_predicates_pos!(verif, nvar, uset)
    for region in uset.regions
        add_predicate!(verif, PosPredicate(nvar, region.domain, region.loc))
    end
end

function _add_predicates_lie!(verif, nvar, sys)
    for piece in sys.pieces
        add_predicate!(verif, LiePredicate(
            nvar, piece.domain, piece.loc1, piece.A, piece.b, piece.loc2
        ))
    end
end

snapshot(::Nothing, ::Any) = nothing

struct TraceRecorder
    mpf_list::Vector{MultiPolyFunc}
    pos_evids_list::Vector{Vector{PosEvidence}}
    lie_evids_list::Vector{Vector{LieEvidence}}
end

TraceRecorder() = TraceRecorder(MultiPolyFunc[], PosEvidence[], LieEvidence[])

function snapshot(tracerec::TraceRecorder, gen::Generator)
    push!(tracerec.pos_evids_list, copy(gen.pos_evids))
    push!(tracerec.lie_evids_list, copy(gen.lie_evids))
end

function snapshot(tracerec::TraceRecorder, mpf::MultiPolyFunc)
    push!(tracerec.mpf_list, mpf)
end

function learn_lyapunov!(
        lear::Learner, iter_max, solver_gen, solver_verif;
        do_print=true, tracerec=nothing
    )
    gen = Generator(lear.nvar, lear.nloc)
    for state in lear.iset.states
        _add_evidences_neg!(gen, state)
    end

    verif = Verifier()
    _add_predicates_pos!(verif, lear.nvar, lear.uset)
    _add_predicates_lie!(verif, lear.nvar, lear.sys)

    mpf = MultiPolyFunc(lear.nloc)
    iter = 0
    params = lear.params
    M, offmax, radmax = params[:bigM], params[:offmax], params[:rmax_gen]
    xmax, objmax = params[:xmax], params[:rmax_verif]
    tol_pos, tol_lie = lear.tols[:pos], lear.tols[:lie]

    while true
        iter += 1
        do_print && println("Iter: ", iter)
        if iter > iter_max
            println(string("Max iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, iter
        end
        snapshot(tracerec, gen)

        # Generator
        mpf, r = compute_mpf_evidence(gen, M, offmax, radmax, solver_gen)
        snapshot(tracerec, mpf)
        if do_print
            println("|--- radius: ", r)
        end
        if r < lear.tols[:rad]
            println(string("Satisfiability radius too small: ", r))
            return RADIUS_TOO_SMALL, mpf, iter
        end

        # Verifier
        do_print && print("|--- Verify pos... ")
        obj, x, loc = verify_pos(verif, mpf, xmax, objmax, solver_verif)
        if obj > tol_pos
            do_print && println("CE found: ", x, ", ", obj, ", ", loc)
            _add_evidences_pos!(gen, State(loc, x))
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify lie... ")
        obj, x, loc = verify_lie(verif, mpf, xmax, objmax, solver_verif)
        if obj > tol_lie
            do_print && println("CE found: ", x, ", ", obj, ", ", loc)
            _add_evidences_lie!(gen, lear.sys, State(loc, x))
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        
        println("No CE found")
        println("Valid CLF: terminated")
        return BARRIER_FOUND, mpf, iter
    end
end