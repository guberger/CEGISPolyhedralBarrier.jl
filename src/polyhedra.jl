struct Halfspace
    a::_VT_
    β::Float64
end

Base.in(x, h::Halfspace) = dot(h.a, x) + h.β ≤ 0

struct Polyhedron
    halfspaces::Vector{Halfspace}
end

Polyhedron() = Polyhedron(Halfspace[])

function add_halfspace!(p::Polyhedron, h::Halfspace)
    push!(p.halfspaces, h)    
end

function add_halfspace!(p::Polyhedron, a, β)
    add_halfspace!(p, Halfspace(a, β))
end

Base.in(x, p::Polyhedron) = all(h -> x ∈ h, p.halfspaces)