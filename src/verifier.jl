abstract type VerifierProblem end

function _optim_verif(
        prob::VerifierProblem,
        afs_dom, afs_safe, afs_inv, afs_BF, af_target,
        A, b, xmax, N, solver
    )
    model = solver()
    x = @variable(model, [1:N], lower_bound=-xmax, upper_bound=xmax)
    r = @variable(model, upper_bound=10*xmax)

    for af in afs_dom
        _add_constr_in_domain!(prob, model, af, x)
    end
    for af in afs_inv
        _add_constr_in_domain!(prob, model, af, x)
    end

    for af in afs_safe
        _add_constr_in_barrier!(prob, model, af, x, r)
    end
    for af in afs_BF
        _add_constr_in_barrier!(prob, model, af, x, r)
    end

    _add_constr_out_target!(prob, model, af_target, A*x + b, r)

    @objective(model, Max, r)

    optimize!(model)

    xopt = has_values(model) ? value.(x) : fill(NaN, N)
    ropt = has_values(model) ? objective_value(model) : -Inf
    ps, ts = primal_status(model), termination_status(model)
    flag = ps == FEASIBLE_POINT && ts == OPTIMAL
    @assert flag || (ps == NO_SOLUTION && ts == INFEASIBLE)

    return xopt, ropt, flag
end

function _verify(
        prob::VerifierProblem,
        mpf_safe, mpf_inv, mpf_BF, mpf_target,
        sys, xmax, N, solver
    )
    xopt::Vector{Float64} = fill(NaN, N)
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
                af_target, A, b, xmax, N, solver
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

function _add_constr_in_domain!(::VerifierProblem, model, af, x)
    @constraint(model, _eval(af, x) ≤ 0)
end

function _add_constr_in_barrier!(::VerifierProblem, model, af, x, r)
    @constraint(model, _eval(af, x) + norm(af.a, Inf)*r ≤ 0)
end

struct VerifierSafe{ET} <: VerifierProblem
    η::ET
end

function _add_constr_out_target!(prob::VerifierSafe, model, af, x, ::Any)
    @constraint(model, _eval(af, x) ≥ norm(af.a, Inf)*prob.η)
end

function verify_safe(
        sys::System, mpf_safe::MultiPolyFunc, mpf_inv::MultiPolyFunc,
        mpf_BF::MultiPolyFunc, xmax, η, N, solver
    )
    prob = VerifierSafe(η)
    return _verify(
        prob, mpf_safe, mpf_inv, mpf_BF, mpf_safe, sys, xmax, N, solver
    )
end

struct VerifierBF{ET} <: VerifierProblem
    η::ET
end

function _add_constr_out_target!(prob::VerifierBF, model, af, x, r)
    @constraint(model, _eval(af, x) ≥ norm(af.a, Inf)*(prob.η + r))
end

function verify_BF(
        sys::System, mpf_safe::MultiPolyFunc, mpf_inv::MultiPolyFunc,
        mpf_BF::MultiPolyFunc, xmax, η, N, solver
    )
    prob = VerifierBF(η)
    return _verify(
        prob, mpf_safe, mpf_inv, mpf_BF, mpf_BF, sys, xmax, N, solver
    )
end