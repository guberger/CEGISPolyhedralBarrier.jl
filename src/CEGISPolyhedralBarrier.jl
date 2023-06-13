module CEGISPolyhedralBarrier

using LinearAlgebra
using JuMP

struct AffForm
    a::Vector{Float64}
    β::Float64
end
_eval(af::AffForm, x) = dot(af.a, x) + af.β
margin(af::AffForm, atol, rtol) = atol + rtol*norm(af.a, Inf)
isless_approx(af::AffForm, x, atol, rtol) = _eval(af, x) ≤ margin(af, atol, rtol)

struct GenForm
    loc::Int
    af::AffForm
end

struct Piece
    afs_dom::Vector{AffForm}
    loc_src::Int
    A::Matrix{Float64}
    b::Vector{Float64}
    loc_dst::Int
end

struct State
    loc::Int
    x::Vector{Float64}
end

struct Link
    src::State
    dst::State
end

include("separator.jl")
include("generator.jl")
include("crosser.jl")
include("verifier.jl")
# include("learner.jl")

end # module
