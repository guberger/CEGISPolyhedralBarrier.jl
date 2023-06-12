module CEGISPolyhedralBarrier

using LinearAlgebra
using JuMP

struct AffForm
    a::Vector{Float64}
    β::Float64
end
_eval(af::AffForm, point) = dot(af.a, point) + af.β
margin(af::AffForm, atol, rtol) = atol + rtol*norm(af.a, Inf)
isless_approx(af::AffForm, point, atol, rtol) =
    _eval(af, point) ≤ margin(af, atol, rtol)

struct PolyFunc
    afs::Vector{AffForm}
end
Base.empty!(pf::PolyFunc) = empty!(pf.afs)
isless_approx(pf::PolyFunc, point, atol, rtol) =
    all(af -> isless_approx(af, point, atol, rtol), pf.afs)
empty_pf() = PolyFunc(AffForm[])

struct MultiPolyFunc
    pfs::Vector{PolyFunc}
end
Base.empty!(mpf::MultiPolyFunc, loc) = empty!(mpf.pfs[loc])
empty_mpf(M) = MultiPolyFunc([empty_pf() for loc = 1:M])

struct Piece
    pf_dom::PolyFunc
    loc1::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    loc2::Int
end

struct System
    pieces::Vector{Piece}
end

struct Grid
    points::Vector{Vector{Float64}}
end
empty_grid() = Grid(Vector{Float64}[])

struct MultiGrid
    grids::Vector{Grid}
end
empty_multigrid(M) = MultiGrid([empty_grid() for loc = 1:M])

include("separator.jl")
include("verifier.jl")
# include("learner.jl")

end # module
