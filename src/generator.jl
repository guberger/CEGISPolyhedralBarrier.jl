struct NegEvidence{N}
    loc::Int
    point::Point{N}
end

struct PosEvidence{N}
    loc::Int
    i::Int
    point::Point{N}
end

struct LieEvidence{N}
    loc1::Int
    i1::Int
    point1::Point{N}
    loc2::Int
    point2::Point{N} # A*point1 + b
    nA::Float64 # opnorm(A)
end

struct Generator{N,M}
    nafs::NTuple{M,Int}
    neg_evids::Vector{NegEvidence{N}}
    pos_evids::Vector{PosEvidence{N}}
    lie_evids::Vector{LieEvidence{N}}
end

Generator{N}(nafs::NTuple{M,Int}) where {N,M} = Generator{N,M}(
    nafs, NegEvidence{N}[], PosEvidence{N}[], LieEvidence{N}[]
)

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

struct _AF{N}
    a::SVector{N,VariableRef}
    β::VariableRef
end
_eval(af::_AF, point) = dot(point, af.a) + af.β

struct _PF{N}
    afs::Vector{_AF{N}}
end
_PF{N}(naf::Int) where N = _PF(Vector{_AF{N}}(undef, naf))

struct _MPF{N,M}
    pfs::NTuple{M,_PF{N}}
end

function _add_vars!(model, ::Val{N}, nafs::NTuple{M,Int}) where {N,M}
    pfs = map(naf -> _PF{N}(naf), nafs)
    for (loc, naf) in enumerate(nafs)
        for i = 1:naf
            a = SVector(ntuple(
                k -> @variable(model, lower_bound=-1, upper_bound=1), Val(N)
            ))
            β = @variable(model, lower_bound=-100, upper_bound=100)
            pfs[loc].afs[i] = _AF(a, β)
        end
    end
    r = @variable(model, upper_bound=10)
    return _MPF(pfs), r
end

function _add_neg_constr!(model, af, r, point, α, γ)
    @constraint(model, _eval(af, point) + α*r + γ ≤ 0)
end

function _add_pos_constr!(model, af, r, bin, point, α, β, γ)
    @constraint(model, _eval(af, point) + β*bin ≥ α*r + γ)
end

function _add_lie_constr(model, af, r, bin, point, α, β, γ)
    @constraint(model, _eval(af, point) + α*r + γ ≤ β*(1 - bin))
end

_value(af::_AF) = AffForm(value.(af.a), value(af.β))

abstract type GeneratorProblem end

function _compute_mpf(
        prob::GeneratorProblem, gen::Generator{N}, solver
    ) where N
    model = solver()
    mpf, r = _add_vars!(model, Val(N), gen.nafs)

    for evid in gen.neg_evids
        for af in mpf.pfs[evid.loc].afs
            _add_constr_prob!(prob, model, af, r, evid)
        end
    end

    for evid in gen.pos_evids
        af = mpf.pfs[evid.loc].afs[evid.i]
        _add_constr_prob!(prob, model, af, r, evid)
    end

    for evid in gen.lie_evids
        bin = @variable(model, binary=true)
        af1 = mpf.pfs[evid.loc1].afs[evid.i1]
        for af2 in mpf.pfs[evid.loc2].afs
            _add_constr_prob!(prob, model, af1, af2, r, bin, evid)
        end
    end

    @objective(model, Max, r)

    optimize!(model)

    @assert termination_status(model) == OPTIMAL
    @assert primal_status(model) == FEASIBLE_POINT

    return MultiPolyFunc(map(
        pf -> PolyFunc([_value(af) for af in pf.afs]), mpf.pfs
    )), value(r)
end

## Feasibility

struct GeneratorFeasibility <: GeneratorProblem
    ϵ::Float64
    M::Float64
end

function _add_constr_prob!(
        prob::GeneratorFeasibility, model, af, r, evid::NegEvidence
    )
    _add_neg_constr!(model, af, r, evid.point, 1, prob.ϵ)
end

function _add_constr_prob!(
        prob::GeneratorFeasibility, model, af, r, evid::PosEvidence
    )
    _add_pos_constr!(model, af, r, 0, evid.point, 1, 0, prob.ϵ)
end

function _add_constr_prob!(
        prob::GeneratorFeasibility, model, af1, af2, r, bin, evid::LieEvidence
    )
    _add_pos_constr!(model, af1, r, bin, evid.point1, 1, prob.M, prob.ϵ)
    _add_lie_constr(model, af2, r, bin, evid.point2, 1, prob.M, prob.ϵ)
end

function compute_mpf_feasibility(gen::Generator, ϵ, M, solver)
    prob = GeneratorFeasibility(ϵ, M)
    return _compute_mpf(prob, gen, solver)
end

## Evidence

struct GeneratorEvidence <: GeneratorProblem
    M::Float64
end

function _add_constr_prob!(
        ::GeneratorEvidence, model, af, r, evid::NegEvidence
    )
    _add_neg_constr!(model, af, r, evid.point, 1, 0)
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

function compute_mpf_evidence(gen::Generator, M, solver)
    prob = GeneratorEvidence(M)
    return _compute_mpf(prob, gen, solver)
end