struct Verifier{N,M}
    mpf_safe::MultiPolyFunc{N,M}
    mpf_inv::MultiPolyFunc{N,M}
    mpf_BF::MultiPolyFunc{N,M}
    sys::System{N}
    xmax::Float64
end

abstract type VerifierProblem end

function _optim_verif(
        prob::VerifierProblem, af1_, af2,
        A::SMatrix{N,N}, b::SVector{N}, xmax, solver
    ) where N
    model = solver()
    x = SVector(ntuple(
        k -> @variable(model, lower_bound=-xmax, upper_bound=xmax), Val(N)
    ))
    r = @variable(model, upper_bound=10*xmax)

    for af1 in af1_
        _add_constr_in!(prob, model, af1, x, r)
    end

    _add_constr_out!(prob, model, af2, A*x + b, r)

    @objective(model, Max, r)

    optimize!(model)

    xopt = has_values(model) ? value.(x) : SVector(ntuple(k -> NaN, Val(N)))
    ropt = has_values(model) ? objective_value(model) : -Inf
    ps, ts = primal_status(model), termination_status(model)
    flag = ps == FEASIBLE_POINT && ts == OPTIMAL
    @assert flag || (ps == NO_SOLUTION && ts == INFEASIBLE)

    return xopt, ropt, flag
end

function _verify(
        prob::VerifierProblem, sys::System{N}, mpfs_in, mpf_out, xmax, solver
    ) where N
    xopt::SVector{N,Float64} = SVector(ntuple(k -> NaN, Val(N)))
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        afs_in_ = map(mpf -> mpf.pfs[piece.loc1].afs, mpfs_in)
        af1_ = Iterators.flatten((pf_dom.afs, afs_in_...))
        for af2 in mpf_out.pfs[piece.loc2].afs
            x, r, flag = _optim_verif(prob, af1_, af2, A, b, xmax, solver)
            if flag && r > ropt
                xopt = x
                ropt = r
                locopt = piece.loc1
            end
        end
    end
    return xopt, ropt, locopt
end

struct VerifierSafe <: VerifierProblem
    δ::Float64
end

function _add_constr_in!(::VerifierSafe, model, af, x, r)
    @constraint(model, _eval(af, x) + norm(af.a, 1)*r ≤ 0)
end

function _add_constr_out!(prob::VerifierSafe, model, af, x, ::Any)
    @constraint(model, _eval(af, x) + prob.δ ≥ 0)
end

function verify_safe(verif::Verifier, δ, solver)
    prob = VerifierSafe(δ)
    mpfs_in = (verif.mpf_safe, verif.mpf_inv, verif.mpf_BF)
    mpf_out = verif.mpf_safe
    return _verify(prob, verif.sys, mpfs_in, mpf_out, verif.xmax, solver)
end

struct VerifierBF <: VerifierProblem
    δ::Float64
end

function _add_constr_in!(::VerifierBF, model, af, x, r)
    @constraint(model, _eval(af, x) + norm(af.a, 1)*r ≤ 0)
end

function _add_constr_out!(prob::VerifierBF, model, af, x, r)
    @constraint(model, _eval(af, x) - norm(af.a, 1)*r + prob.δ ≥ 0)
end

function verify_BF(verif::Verifier, δ, solver)
    prob = VerifierBF(δ)
    mpfs_in = (verif.mpf_safe, verif.mpf_inv, verif.mpf_BF)
    mpf_out = verif.mpf_BF
    return _verify(prob, verif.sys, mpfs_in, mpf_out, verif.xmax, solver)
end