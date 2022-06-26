module CEGISPolyhedralBarrier

using LinearAlgebra
using JuMP

_RSC_ = JuMP.MathOptInterface.ResultStatusCode
_TSC_ = JuMP.MathOptInterface.TerminationStatusCode
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

struct MultiPolyFunc
    pfs::Vector{PolyFunc}
end

MultiPolyFunc(nloc::Int) = MultiPolyFunc([PolyFunc() for loc = 1:nloc])

# struct Piece
#     domain::Polyhedron
#     loc1::Int
#     A::_MT_
#     loc2::Int
# end

# struct System
#     pieces::Vector{Piece}
# end

# System() = System(Piece[])

# function add_piece!(sys::System, domain, loc1, A, loc2)
#     push!(sys.pieces, Piece(domain, loc1, A, loc2))
# end

include("generator.jl")
# include("verifier.jl")
# include("learner.jl")

end # module
