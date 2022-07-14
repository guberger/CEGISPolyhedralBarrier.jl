struct Generator{N,M}
    neg_evids::PointSet{N,M}
    pos_evids::PointSet{N,M}
end

Generator{N,M}() where {N,M} = Generator(PointSet{N,M}(), PointSet{N,M}())

function add_neg_evidence!(gen::Generator, loc, point)
    add_point!(gen.neg_evids, loc, point)
end

function add_pos_evidence!(gen::Generator, loc, point)
    add_point!(gen.pos_evids, loc, point)
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

function _add_vars!(model, ::Val{N}, βmax, nafs::NTuple{M,Int}) where {N,M}
    pfs = map(naf -> _PF{N}(naf), nafs)
    for (loc, naf) in enumerate(nafs)
        for i = 1:naf
            a = SVector(ntuple(
                k -> @variable(model, lower_bound=-1, upper_bound=1), Val(N)
            ))
            β = @variable(model, lower_bound=-βmax, upper_bound=βmax)
            pfs[loc].afs[i] = _AF(a, β)
        end
    end
    r = @variable(model, upper_bound=10)
    return _MPF(pfs), r
end

function _add_neg_constr!(model, af, r, point, α, γ)
    @constraint(model, _eval(af, point) + α*r + γ ≤ 0)
end

function _add_pos_constr!(model, af, r, point, α, γ)
    @constraint(model, _eval(af, point) - α*r - γ ≥ 0)
end

_value(af::_AF) = AffForm(value.(af.a), value(af.β))

abstract type GeneratorProblem end

function _compute_mpf(
        prob::GeneratorProblem, gen::Generator{N}, βmax, solver
    ) where N
    model = solver()
    nafs = length.(gen.pos_evids.points_list)
    mpf, r = _add_vars!(model, Val(N), βmax, nafs)

    for (pf, points) in zip(mpf.pfs, gen.neg_evids.points_list)
        for (af, point) in Iterators.product(pf.afs, points)
            _add_constr_neg_evid!(prob, model, af, r, point)
        end
    end

    for (pf, points) in zip(mpf.pfs, gen.pos_evids.points_list)
        for (af, point) in zip(pf.afs, points)
            _add_constr_pos_evid!(prob, model, af, r, point)
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

## Robust

struct GeneratorRobust <: GeneratorProblem
    ϵ::Float64
end

function _add_constr_neg_evid!(prob::GeneratorRobust, model, af, r, point)
    _add_neg_constr!(model, af, r, point, 1, prob.ϵ)
end

function _add_constr_pos_evid!(prob::GeneratorRobust, model, af, r, point)
    _add_pos_constr!(model, af, r, point, 1, prob.ϵ)
end

function compute_mpf_robust(gen::Generator, ϵ, βmax, solver)
    prob = GeneratorRobust(ϵ)
    return _compute_mpf(prob, gen, βmax, solver)
end