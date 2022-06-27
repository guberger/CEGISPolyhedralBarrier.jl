module CEGISPolyhedralBarrier

using LinearAlgebra
using JuMP

_status(model) = (primal_status(model), termination_status(model))
_VT_ = Vector{Float64}
_MT_ = Matrix{Float64}
Point = _VT_

include("polyhedra.jl")

struct AffForm
    lin::_VT_
    off::Float64
end
_eval(af::AffForm, point) = dot(af.lin, point) + af.off

struct PolyFunc
    afs::Vector{AffForm}
end

PolyFunc() = PolyFunc(AffForm[])
add_af!(pf::PolyFunc, af::AffForm) = push!(pf.afs, af)

struct MultiPolyFunc
    pfs::Vector{PolyFunc}
end

MultiPolyFunc(nloc::Int) = MultiPolyFunc([PolyFunc() for loc = 1:nloc])
add_af!(mpf::MultiPolyFunc, loc::Int, af::AffForm) = add_af!(mpf.pfs[loc], af)

struct Piece
    domain::Polyhedron
    loc1::Int
    A::_MT_
    b::_VT_
    loc2::Int
end

struct System
    pieces::Vector{Piece}
end

System() = System(Piece[])

function add_piece!(sys::System, domain, loc1, A, b, loc2)
    push!(sys.pieces, Piece(domain, loc1, A, b, loc2))
end

struct State
    loc::Int
    point::Point
end

struct Region
    loc::Int
    domain::Polyhedron
end

struct InitSet
    states::Vector{State}
end

InitSet() = InitSet(State[])
add_state!(iset::InitSet, state::State) = push!(iset.states, state)

struct UnsafeSet
    regions::Vector{Region}
end

UnsafeSet() = UnsafeSet(Region[])
add_region!(uset::UnsafeSet, region::Region) = push!(uset.regions, region)

include("generator.jl")
include("verifier.jl")
include("learner.jl")

end # module
