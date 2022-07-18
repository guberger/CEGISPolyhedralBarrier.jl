abstract type VerifierProblem end

function _optim_verif(
        prob::VerifierProblem, afs_dom, afs_safe, afs_inv, afs_BF,
        af_target, A::SMatrix{N,N}, b::SVector{N}, xmax, solver
    ) where N
    model = solver()
    x = SVector(ntuple(
        k -> @variable(model, lower_bound=-xmax, upper_bound=xmax), Val(N)
    ))
    r = @variable(model, upper_bound=10*xmax)

    for af in Iterators.flatten((afs_dom, afs_inv))
        _add_constr_in_domain!(prob, model, af, x)
    end

    for af in Iterators.flatten((afs_safe, afs_BF))
        _add_constr_in_barrier!(prob, model, af, x, r)
    end

    _add_constr_out_target!(prob, model, af_target, A*x + b, r)

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
        prob::VerifierProblem, sys::System{N},
        mpf_safe, mpf_inv, mpf_BF,
        mpf_target, xmax, solver
    ) where N
    xopt::SVector{N,Float64} = SVector(ntuple(k -> NaN, Val(N)))
    ropt::Float64 = -Inf
    locopt::Int = 0
    for piece in sys.pieces
        pf_dom, A, b = piece.pf_dom, piece.A, piece.b
        afs_dom = pf_dom.afs
        afs_safe = mpf_safe.pfs[piece.loc1].afs
        afs_inv = mpf_inv.pfs[piece.loc1].afs
        afs_BF = mpf_BF.pfs[piece.loc1].afs
        for af_target in mpf_target.pfs[piece.loc2].afs
            x, r, flag = _optim_verif(
                prob, afs_dom, afs_safe, afs_inv, afs_BF,
                af_target, A, b, xmax, solver
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

struct VerifierSafe <: VerifierProblem
    η::Float64
end

function _add_constr_in_domain!(::VerifierSafe, model, af, x)
    @constraint(model, _eval(af, x) ≤ 0)
end

function _add_constr_in_barrier!(::VerifierSafe, model, af, x, r)
    @constraint(model, _eval(af, x) + norm(af.a, Inf)*r ≤ 0)
end

function _add_constr_out_target!(prob::VerifierSafe, model, af, x, ::Any)
    @constraint(model, _eval(af, x) ≥ norm(af.a, Inf)*prob.η)
end

function verify_safe(
        sys::System, mpf_safe::MultiPolyFunc, mpf_inv::MultiPolyFunc,
        mpf_BF::MultiPolyFunc, xmax::Float64, η, solver
    )
    prob = VerifierSafe(η)
    return _verify(
        prob, sys, mpf_safe, mpf_inv, mpf_BF, mpf_safe, xmax, solver
    )
end

struct VerifierBF <: VerifierProblem
    η::Float64
end

function _add_constr_in_domain!(::VerifierBF, model, af, x)
    @constraint(model, _eval(af, x) ≤ 0)
end

function _add_constr_in_barrier!(::VerifierBF, model, af, x, r)
    @constraint(model, _eval(af, x) + norm(af.a, Inf)*r ≤ 0)
end

function _add_constr_out_target!(prob::VerifierBF, model, af, x, r)
    @constraint(model, _eval(af, x) ≥ norm(af.a, Inf)*(prob.η + r))
end

function verify_BF(
        sys::System, mpf_safe::MultiPolyFunc, mpf_inv::MultiPolyFunc,
        mpf_BF::MultiPolyFunc, xmax::Float64, η, solver
    )
    prob = VerifierBF(η)
    return _verify(
        prob, sys, mpf_safe, mpf_inv, mpf_BF, mpf_BF, xmax, solver
    )
end