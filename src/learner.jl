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
    tols::Dict{Symbol,Float64}
    params::Dict{Symbol,Float64}
end

function Learner(sys, mpf_safe, mpf_inv, isetϵ, ϵ)
    tols = Dict([
        :rad => eps(1.0),
        :verif => -eps(1.0),
        :dom => 1e-8
    ])
    params = Dict([
        :βmax => 1e3,
        :bigM => 1e3,
        :xmax => 1e3,
    ])
    return Learner(sys, mpf_safe, mpf_inv, isetϵ, ϵ, tols, params)
end

_setsafe!(D, k, v) = (@assert haskey(D, k); D[k] = v)
set_tol!(lear::Learner, s::Symbol, v) = _setsafe!(lear.tols, s, v)
set_param!(lear::Learner, s::Symbol, v) = _setsafe!(lear.params, s, v)

## Learn Barrier
function _add_evidences_pos(gen, loc, point)
    nafs = setindex(gen.nafs, gen.nafs[loc] + 1, loc)
    gen = Generator(nafs, gen.neg_evids, gen.pos_evids, gen.lie_evids)
    add_evidence!(gen, PosEvidence(loc, nafs[loc], point))
    return gen
end

function _add_evidences_lie(gen, sys, loc, point, tol_dom)
    nafs = setindex(gen.nafs, gen.nafs[loc] + 1, loc)
    gen = Generator(nafs, gen.neg_evids, gen.pos_evids, gen.lie_evids)
    i = nafs[loc]
    for piece in sys.pieces
        loc != piece.loc1 && continue
        !_neg(piece.pf_dom, point, tol_dom) && continue
        point2 = piece.A*point + piece.b
        loc2 = piece.loc2
        nA = opnorm(piece.A, 1)
        add_evidence!(gen, LieEvidence(loc, i, point, loc2, point2, nA))
    end
    return gen
end

function learn_lyapunov!(
        lear::Learner{N,M}, iter_max, solver_gen, solver_verif;
        do_print=true, tracerec=nothing
    ) where {N,M}
    gen = Generator{N}(ntuple(loc -> 0, Val(M)))
    for (loc, points) in enumerate(lear.iset.points_list)
        for point in points
            add_evidence!(gen, NegEvidence(loc, point))
        end
    end

    mpf = MultiPolyFunc{N,M}()
    iter = 0
    params = lear.params
    βmax, bigM, xmax = params[:βmax], params[:bigM], params[:xmax]
    tol_dom = lear.tols[:dom]

    while true
        iter += 1
        do_print && println("Iter: ", iter)
        if iter > iter_max
            println(string("Max iter exceeded: ", iter))
            return MAX_ITER_REACHED, mpf, gen, iter
        end

        # Generator
        mpf, r = compute_mpf_evidence(gen, bigM, βmax, solver_gen)
        if do_print
            println("|--- radius: ", r)
        end
        if r < lear.tols[:rad]
            println(string("Satisfiability radius too small: ", r))
            return RADIUS_TOO_SMALL, mpf, gen, iter
        end

        # Verifier
        verif = Verifier(lear.mpf_safe, lear.mpf_inv, mpf, lear.sys, xmax)
        do_print && print("|--- Verify safe... ")
        x, obj, loc = verify_safe(verif, solver_verif)
        if obj > lear.tols[:verif]
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            gen = _add_evidences_pos(gen, loc, x)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        do_print && print("|--- Verify BF... ")
        x, obj, loc = verify_BF(verif, solver_verif)
        if obj > lear.tols[:verif]
            do_print && println("CE found: ", x, ", ", loc, ", ", obj)
            gen = _add_evidences_lie(gen, lear.sys, loc, x, tol_dom)
            continue
        else
            do_print && println("No CE found: ", obj)
        end
        
        println("Valid CLF: terminated")
        return BARRIER_FOUND, mpf, gen, iter
    end
end