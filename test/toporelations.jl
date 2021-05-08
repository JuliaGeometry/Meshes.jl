@testset "TopologicalRelation" begin
  @testset "HalfEdgeStructure" begin
    # 2 triangles
    elems = connect.([(1,2,3),(4,3,2)])
    struc = HalfEdgeStructure(elems)
    ∂ = Boundary{2,0}(struc)
    @test ∂(1) == [2,3,1]
    @test ∂(2) == [3,2,4]
    ∂ = Boundary{2,1}(struc)
    @test ∂(1) == connect.([(2,3),(3,1),(1,2)])
    @test ∂(2) == connect.([(3,2),(2,4),(4,3)])
    𝒞 = Coboundary{0,1}(struc)
    @test 𝒞(1) == connect.([(1,2),(1,3)])
    @test 𝒞(2) == connect.([(2,4),(2,3),(2,1)])
    @test 𝒞(3) == connect.([(3,1),(3,2),(3,4)])
    @test 𝒞(4) == connect.([(4,3),(4,2)])
    𝒞 = Coboundary{0,2}(struc)
    @test 𝒞(1) == connect.([(1,2,3)])
    @test 𝒞(2) == connect.([(2,4,3),(2,3,1)])
    @test 𝒞(3) == connect.([(3,1,2),(3,2,4)])
    @test 𝒞(4) == connect.([(4,3,2)])
    𝒞 = Coboundary{1,2}(struc)
    @test 𝒞(1) == connect.([(3,2,4),(2,3,1)])
    @test 𝒞(2) == connect.([(1,2,3)])
    @test 𝒞(3) == connect.([(3,1,2)])
    @test 𝒞(4) == connect.([(2,4,3)])
    @test 𝒞(5) == connect.([(4,3,2)])
    𝒜 = Adjacency{0}(struc)
    @test 𝒜(1) == [2,3]
    @test 𝒜(2) == [4,3,1]
    @test 𝒜(3) == [1,2,4]
    @test 𝒜(4) == [3,2]

    # 2 triangles + 2 quadrangles
    elems = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)])
    struc = HalfEdgeStructure(elems)
    ∂ = Boundary{2,0}(struc)
    @test ∂(1) == [1,2,6,5]
    @test ∂(2) == [6,2,4]
    @test ∂(3) == [6,4,3,5]
    @test ∂(4) == [3,1,5]
    ∂ = Boundary{2,1}(struc)
    @test ∂(1) == connect.([(1,2),(2,6),(6,5),(5,1)])
    @test ∂(2) == connect.([(6,2),(2,4),(4,6)])
    @test ∂(3) == connect.([(6,4),(4,3),(3,5),(5,6)])
    @test ∂(4) == connect.([(3,1),(1,5),(5,3)])
    𝒞 = Coboundary{0,1}(struc)
    @test 𝒞(1) == connect.([(1,2),(1,5),(1,3)])
    @test 𝒞(2) == connect.([(2,4),(2,6),(2,1)])
    @test 𝒞(3) == connect.([(3,1),(3,5),(3,4)])
    @test 𝒞(4) == connect.([(4,3),(4,6),(4,2)])
    @test 𝒞(5) == connect.([(5,6),(5,3),(5,1)])
    @test 𝒞(6) == connect.([(6,2),(6,4),(6,5)])
    𝒞 = Coboundary{0,2}(struc)
    @test 𝒞(1) == connect.([(1,2,6,5),(1,5,3)])
    @test 𝒞(2) == connect.([(2,4,6),(2,6,5,1)])
    @test 𝒞(3) == connect.([(3,1,5),(3,5,6,4)])
    @test 𝒞(4) == connect.([(4,3,5,6),(4,6,2)])
    @test 𝒞(5) == connect.([(5,6,4,3),(6,5,1,2),(5,3,1)])
    @test 𝒞(6) == connect.([(6,2,4),(2,6,5,1),(6,4,3,5)])
    𝒞 = Coboundary{1,2}(struc)
    @test 𝒞(1) == connect.([(1,2,6,5)])
    @test 𝒞(2) == connect.([(3,1,5)])
    @test 𝒞(3) == connect.([(6,2,4),(2,6,5,1)])
    @test 𝒞(4) == connect.([(4,6,2),(6,4,3,5)])
    @test 𝒞(5) == connect.([(5,6,4,3),(6,5,1,2)])
    @test 𝒞(6) == connect.([(1,5,3),(5,1,2,6)])
    @test 𝒞(7) == connect.([(4,3,5,6)])
    @test 𝒞(8) == connect.([(3,5,6,4),(5,3,1)])
    @test 𝒞(9) == connect.([(2,4,6)])
    𝒜 = Adjacency{0}(struc)
    @test 𝒜(1) == [2,5,3]
    @test 𝒜(2) == [4,6,1]
    @test 𝒜(3) == [1,5,4]
    @test 𝒜(4) == [3,6,2]
    @test 𝒜(5) == [6,3,1]
    @test 𝒜(6) == [2,4,5]
  end
end
