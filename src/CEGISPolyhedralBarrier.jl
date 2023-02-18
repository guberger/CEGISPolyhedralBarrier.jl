module CEGISPolyhedralBarrier

using LinearAlgebra
using StaticArrays
using JuMP

Point{N} = SVector{N,Float64}

struct AffForm{N}
    a::Point{N}
    β::Float64
end
_eval(af::AffForm, point) = dot(af.a, point) + af.β

struct PolyFunc{N}
    afs::Vector{AffForm{N}}
end

PolyFunc{N}() where N = PolyFunc(AffForm{N}[])
add_af!(pf::PolyFunc, af::AffForm) = push!(pf.afs, af)
add_af!(pf::PolyFunc, af_...) = add_af!(pf, AffForm(af_...))
Base.empty!(pf::PolyFunc) = empty!(pf.afs)
_prox(pf::PolyFunc, point, r) = all(
    af -> _eval(af, point) ≤ norm(af.a, Inf)*r, pf.afs
)

struct MultiPolyFunc{N,M}
    pfs::NTuple{M,PolyFunc{N}}
end

MultiPolyFunc{N,M}() where {N,M} = MultiPolyFunc(
    ntuple(loc -> PolyFunc{N}(), Val(M))
)
add_af!(mpf::MultiPolyFunc, loc::Int, af_...) = add_af!(mpf.pfs[loc], af_...)
Base.empty!(mpf::MultiPolyFunc, loc::Int) = empty!(mpf.pfs[loc])

struct Piece{N}
    pf_dom::PolyFunc{N}
    loc1::Int
    A::SMatrix{N,N,Float64}
    b::SVector{N,Float64}
    loc2::Int
end

struct System{N}
    pieces::Vector{Piece{N}}
end

System{N}() where N = System(Piece{N}[])
add_piece!(sys::System, piece::Piece) = push!(sys.pieces, piece)
add_piece!(sys::System, piece_...) = add_piece!(sys, Piece(piece_...))

struct MultiSet{N,M}
    sets::NTuple{M,Vector{Point{N}}}
end

MultiSet{N,M}() where {N,M} = MultiSet(ntuple(loc -> Point{N}[], Val(M)))
add_point!(S::MultiSet, loc::Int, point) = push!(S.sets[loc], point)
pop_point!(S::MultiSet, loc::Int) = pop!(S.sets[loc])
Base.empty!(S::MultiSet, loc::Int) = empty!(S.sets[loc])

include("separator.jl")
include("verifier.jl")
include("learner.jl")

end # module
