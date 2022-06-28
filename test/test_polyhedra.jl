using Test
@static if isdefined(Main, :TestLocal)
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier

p1 = CPB.Polyhedron()
CPB.add_halfspace!(p1, [1, 1], 1)

@testset "polyhedron" begin
    @test [1, -2.1] ∈ p1
end

p2 = CPB.Polyhedron()
CPB.add_halfspace!(p2, [-1, -1], -2)
p3 = p1 ∩ p2

@testset "polyhedron" begin
    @test [0, -2.1] ∉ p3
end

nothing