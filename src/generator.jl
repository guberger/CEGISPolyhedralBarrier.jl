struct PosEvidence
    loc::Int
    i::Int
    point::Point
end

struct NegEvidence
    loc::Int
    point::Point
end

struct LieEvidence
    loc1::Int
    i1::Int
    point1::Point
    loc2::Int
    point2::Point # A*point1 + b
    nA::Float64 # opnorm(A)
end

struct Generator
    nvar::Int
    nafs::Vector{Int}
    neg_evids::Vector{NegEvidence}
    pos_evids::Vector{PosEvidence}
    lie_evids::Vector{LieEvidence}
end

Generator(nvar::Int, nloc::Int) = Generator(
    nvar, zeros(Int, nloc), NegEvidence[], PosEvidence[], LieEvidence[]
)

function add_af!(gen::Generator, loc)
    return gen.nafs[loc] += 1
end

function add_evidence!(gen::Generator, evid::NegEvidence)
    push!(gen.neg_evids, evid)
end

function add_evidence!(gen::Generator, evid::PosEvidence)
    push!(gen.pos_evids, evid)
end

function add_evidence!(gen::Generator, evid::LieEvidence)
    push!(gen.lie_evids, evid)
end

## Compute afs

struct _AF
    lin::Vector{VariableRef}
    off::VariableRef
end
_eval(af::_AF, point) = dot(point, af.lin) + af.off

struct _PF
    afs::Vector{_AF}
end

function _add_vars!(model, nvar, nafs, offmax, rmax)
    pfs = Vector{_PF}(undef, length(nafs))
    for (loc, naf) in enumerate(nafs)
        pfs[loc] = _PF(Vector{_AF}(undef, naf))
        for i = 1:naf
            lin = @variable(model, [1:nvar], lower_bound=-1, upper_bound=1)
            off = @variable(model, lower_bound=-offmax, upper_bound=offmax)
            pfs[loc].afs[i] = _AF(lin, off)
            alin = @variable(model, [1:nvar], lower_bound=0, upper_bound=1)
            @constraint(model, -lin .≤ alin)
            @constraint(model, +lin .≤ alin)
            @constraint(model, sum(alin) ≤ 1)
        end
    end
    r = @variable(model, upper_bound=rmax)
    return pfs, r
end

function _add_neg_constr!(model, af, point)
    @constraint(model, _eval(af, point) ≤ 0)
end

function _add_pos_constr!(model, af, r, bin, point, α, β, γ)
    @constraint(model, _eval(af, point) + β*bin ≥ α*r + γ)
end

function _add_lie_constr(model, af, r, bin, point, α, β, γ)
    @constraint(model, _eval(af, point) + α*r + γ ≤ β*(1 - bin))
end

_value(af::_AF) = AffForm(value.(af.lin), value(af.off))

abstract type GeneratorProblem end

_compute_mpf(
    prob::GeneratorProblem, gen::Generator, offmax, rmax, solver
) = _compute_mpf(
    prob, gen.nvar, gen.nafs,
    gen.neg_evids, gen.pos_evids, gen.lie_evids,
    offmax, rmax, solver
)

function _compute_mpf(
        prob::GeneratorProblem, nvar, nafs,
        neg_evids, pos_evids, lie_evids,
        offmax, rmax, solver
    )
    model = solver()
    pfs, r = _add_vars!(model, nvar, nafs, offmax, rmax)

    for evid in neg_evids
        for af in pfs[evid.loc].afs
            _add_constr_prob!(prob, model, af, evid)
        end
    end

    bins = []

    for evid in pos_evids
        af = pfs[evid.loc].afs[evid.i]
        _add_constr_prob!(prob, model, af, r, evid)
    end

    for evid in lie_evids
        bin = @variable(model, binary=true)
        push!(bins, bin)
        af1 = pfs[evid.loc1].afs[evid.i1]
        for af2 in pfs[evid.loc2].afs
            _add_constr_prob!(prob, model, af1, af2, r, bin, evid)
        end
    end

    @objective(model, Max, r)

    optimize!(model)

    if !(_status(model) == (FEASIBLE_POINT, OPTIMAL))
        error(string(
            "Generator: not optimal: ",
            primal_status(model), " ",
            dual_status(model), " ",
            termination_status(model)
        ))
    end

    return MultiPolyFunc([
        PolyFunc([_value(af) for af in pf.afs]) for pf in pfs
    ]), value(r)
end

function _add_constr_prob!(
        ::GeneratorProblem, model, af, evid::NegEvidence
    )
    _add_neg_constr!(model, af, evid.point)
end

## Feasibility

struct GeneratorFeasibility <: GeneratorProblem
    ϵ::Float64
    δ::Float64
    M::Float64
end

function _add_constr_prob!(
        prob::GeneratorFeasibility, model, af, r, evid::PosEvidence
    )
    _add_pos_constr!(model, af, r, 0, evid.point, 1, 0, 1/prob.ϵ)
end

function _add_constr_prob!(
        prob::GeneratorFeasibility, model, af1, af2, r, bin, evid::LieEvidence
    )
    _add_pos_constr!(model, af1, r, bin, evid.point1, 1, prob.M, prob.δ)
    _add_lie_constr(model, af2, r, bin, evid.point2, 1, prob.M, prob.δ)
end

function compute_mpf_feasibility(gen::Generator, ϵ, δ, M, offmax, rmax, solver)
    prob = GeneratorFeasibility(ϵ, δ, M)
    return _compute_mpf(prob, gen, offmax, rmax, solver)
end

## Evidence

struct GeneratorEvidence <: GeneratorProblem
    M::Float64
end

function _add_constr_prob!(
        ::GeneratorEvidence, model, af, r, evid::PosEvidence
    )
    _add_pos_constr!(model, af, r, 0, evid.point, 1, 0, 0)
end

function _add_constr_prob!(
        prob::GeneratorEvidence, model, af1, af2, r, bin, evid::LieEvidence
    )
    _add_pos_constr!(model, af1, r, bin, evid.point1, 1, prob.M, 0)
    _add_lie_constr(model, af2, r, bin, evid.point2, evid.nA, prob.M, 0)
end

function compute_mpf_evidence(gen::Generator, M, offmax, rmax, solver)
    prob = GeneratorEvidence(M)
    return _compute_mpf(prob, gen, offmax, rmax, solver)
end