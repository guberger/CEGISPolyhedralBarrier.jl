using Test
@static if isdefined(Main, :TestLocal)
    include("../src/CEGISPolyhedralBarrier.jl")
else
    using CEGISPolyhedralBarrier
end
CPB = CEGISPolyhedralBarrier

poly = CPB.Polyhedron()
CPB.add_halfspace!(poly, [1.0, 1.0], 1.0)

@testset "polyhedron" begin
    @test [1, -2.1] ∈ poly
end

nothing