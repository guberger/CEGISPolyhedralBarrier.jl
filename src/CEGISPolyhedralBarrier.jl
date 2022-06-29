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
add_af!(pf::PolyFunc, lin, off) = add_af!(pf, AffForm(lin, off))

struct MultiPolyFunc
    pfs::Vector{PolyFunc}
end

MultiPolyFunc(nloc::Int) = MultiPolyFunc([PolyFunc() for loc = 1:nloc])
add_af!(mpf::MultiPolyFunc, loc::Int, af_...) = add_af!(mpf.pfs[loc], af_...)

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
add_piece!(sys::System, piece::Piece) = push!(sys.pieces, piece)
add_piece!(sys::System, domain, loc1, A, b, loc2) = add_piece!(
    sys, Piece(domain, loc1, A, b, loc2)
)

struct State
    loc::Int
    point::Point
end

struct InitialSet
    states::Vector{State}
end

InitialSet() = InitialSet(State[])
add_state!(iset::InitialSet, state::State) = push!(iset.states, state)
add_state!(iset::InitialSet, loc, point) = add_state!(iset, State(loc, point))

struct Region
    loc::Int
    domain::Polyhedron
end

struct UnsafeSet
    regions::Vector{Region}
end

UnsafeSet() = UnsafeSet(Region[])
add_region!(uset::UnsafeSet, region::Region) = push!(uset.regions, region)
add_region!(uset::UnsafeSet, loc, domain) = add_region!(
    uset, Region(loc, domain)
)

include("generator.jl")
include("verifier.jl")
include("learner.jl")

end # module
