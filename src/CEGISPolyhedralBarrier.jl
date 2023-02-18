module CEGISPolyhedralBarrier

using LinearAlgebra
using JuMP

struct AffForm{AT<:AbstractVector,BT}
    a::AT
    β::BT
end
_eval(af::AffForm, point) = dot(af.a, point) + af.β

struct PolyFunc{AT<:AffForm}
    afs::Vector{AT}
end

Base.empty!(pf::PolyFunc) = empty!(pf.afs)
_prox(pf::PolyFunc, point, r) = all(
    af -> _eval(af, point) ≤ norm(af.a, Inf)*r, pf.afs
)

struct MultiPolyFunc{PT<:PolyFunc}
    pfs::Vector{PT}
end

Base.empty!(mpf::MultiPolyFunc, loc) = empty!(mpf.pfs[loc])

struct Piece{PT<:PolyFunc,AT<:AbstractMatrix,BT<:AbstractVector}
    pf_dom::PT
    loc1::Int
    A::AT
    b::BT
    loc2::Int
end

struct System{PT<:Piece}
    pieces::Vector{PT}
end

include("separator.jl")
include("verifier.jl")
include("learner.jl")

end # module
