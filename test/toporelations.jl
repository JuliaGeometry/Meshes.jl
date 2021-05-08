@testset "TopologicalRelation" begin
  @testset "HalfEdgeStructure" begin
    # 2 triangles
    elems = connect.([(1,2,3),(4,3,2)])
    struc = HalfEdgeStructure(elems)
    ∂ₒ = Boundary{2,0}(struc)
    @test ∂ₒ(1) == [2,3,1]
    @test ∂ₒ(2) == [3,2,4]
    ∂₁ = Boundary{2,1}(struc)
    @test ∂₁(1) == connect.([(2,3),(3,1),(1,2)])
    @test ∂₁(2) == connect.([(3,2),(2,4),(4,3)])
    𝒞₁ = Coboundary{0,1}(struc)
    @test 𝒞₁(1) == connect.([(1,2),(1,3)])
    @test 𝒞₁(2) == connect.([(2,4),(2,3),(2,1)])
    @test 𝒞₁(3) == connect.([(3,1),(3,2),(3,4)])
    @test 𝒞₁(4) == connect.([(4,3),(4,2)])
    𝒞₂ = Coboundary{0,2}(struc)
    @test 𝒞₂(1) == connect.([(1,2,3)])
    @test 𝒞₂(2) == connect.([(2,4,3),(2,3,1)])
    @test 𝒞₂(3) == connect.([(3,1,2),(3,2,4)])
    @test 𝒞₂(4) == connect.([(4,3,2)])
  end
end
