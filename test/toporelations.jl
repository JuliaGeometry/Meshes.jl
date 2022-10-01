@testset "TopologicalRelation" begin
  @testset "GridTopology" begin
    # 3 segments
    t = GridTopology(3)
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1,2]
    @test ∂(2) == [2,3]
    @test ∂(3) == [3,4]

    # 2x3 quadrangles
    t = GridTopology(2, 3)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1,2,5,4]
    @test ∂(2) == [2,3,6,5]
    @test ∂(3) == [4,5,8,7]
    @test ∂(4) == [5,6,9,8]
    @test ∂(5) == [7,8,11,10]
    @test ∂(6) == [8,9,12,11]
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1,2]
    @test ∂(2) == [2,3]
    @test ∂(3) == [4,5]
    @test ∂(4) == [5,6]
    @test ∂(5) == [7,8]
    @test ∂(6) == [8,9]
    @test ∂(7) == [10,11]
    @test ∂(8) == [11,12]
    @test ∂(9) == [1,4]
    @test ∂(10) == [2,5]
    @test ∂(11) == [3,6]
    @test ∂(12) == [4,7]
    @test ∂(13) == [5,8]
    @test ∂(14) == [6,9]
    @test ∂(15) == [7,10]
    @test ∂(16) == [8,11]
    @test ∂(17) == [9,12]

    # 2x3x2 hexahedrons
    t = GridTopology(2, 3, 2)
    ∂ = Boundary{3,0}(t)
    @test ∂(1) == [1,2,5,4,13,14,17,16]
    @test ∂(2) == [2,3,6,5,14,15,18,17]
    @test ∂(3) == [4,5,8,7,16,17,20,19]
    @test ∂(12) == [20,21,24,23,32,33,36,35]

    # invalid relations
    t = GridTopology(2, 3)
    @test_throws AssertionError Boundary{3,0}(t)
    @test_throws AssertionError Coboundary{0,3}(t)
    @test_throws AssertionError Adjacency{3}(t)
    @test_throws AssertionError Boundary{0,2}(t)
    @test_throws AssertionError Coboundary{2,0}(t)
  end

  @testset "HalfEdgeTopology" begin
    # 2 triangles
    elems = connect.([(1,2,3),(4,3,2)])
    t = HalfEdgeTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [2,3,1]
    @test ∂(2) == [3,2,4]
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1,3,2]
    @test ∂(2) == [1,4,5]
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [3,2]
    @test ∂(2) == [1,2]
    @test ∂(3) == [3,1]
    @test ∂(4) == [2,4]
    @test ∂(5) == [4,3]
    𝒞 = Coboundary{0,1}(t)
    @test 𝒞(1) == [2,3]
    @test 𝒞(2) == [4,1,2]
    @test 𝒞(3) == [3,1,5]
    @test 𝒞(4) == [5,4]
    𝒞 = Coboundary{0,2}(t)
    @test 𝒞(1) == [1]
    @test 𝒞(2) == [2,1]
    @test 𝒞(3) == [1,2]
    @test 𝒞(4) == [2]
    𝒞 = Coboundary{1,2}(t)
    @test 𝒞(1) == [2,1]
    @test 𝒞(2) == [1]
    @test 𝒞(3) == [1]
    @test 𝒞(4) == [2]
    @test 𝒞(5) == [2]
    𝒜 = Adjacency{0}(t)
    @test 𝒜(1) == [2,3]
    @test 𝒜(2) == [4,3,1]
    @test 𝒜(3) == [1,2,4]
    @test 𝒜(4) == [3,2]

    # 2 triangles + 2 quadrangles
    elems = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(1,5,3)])
    t = HalfEdgeTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1,2,6,5]
    @test ∂(2) == [6,2,4]
    @test ∂(3) == [6,4,3,5]
    @test ∂(4) == [3,1,5]
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1,3,5,6]
    @test ∂(2) == [3,9,4]
    @test ∂(3) == [4,7,8,5]
    @test ∂(4) == [2,6,8]
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1,2]
    @test ∂(2) == [3,1]
    @test ∂(3) == [6,2]
    @test ∂(4) == [4,6]
    @test ∂(5) == [5,6]
    @test ∂(6) == [1,5]
    @test ∂(7) == [4,3]
    @test ∂(8) == [3,5]
    @test ∂(9) == [2,4]
    𝒞 = Coboundary{0,1}(t)
    @test 𝒞(1) == [1,6,2]
    @test 𝒞(2) == [9,3,1]
    @test 𝒞(3) == [2,8,7]
    @test 𝒞(4) == [7,4,9]
    @test 𝒞(5) == [5,8,6]
    @test 𝒞(6) == [3,4,5]
    𝒞 = Coboundary{0,2}(t)
    @test 𝒞(1) == [1,4]
    @test 𝒞(2) == [2,1]
    @test 𝒞(3) == [4,3]
    @test 𝒞(4) == [3,2]
    @test 𝒞(5) == [3,4,1]
    @test 𝒞(6) == [2,3,1]
    𝒞 = Coboundary{1,2}(t)
    @test 𝒞(1) == [1]
    @test 𝒞(2) == [4]
    @test 𝒞(3) == [2,1]
    @test 𝒞(4) == [2,3]
    @test 𝒞(5) == [3,1]
    @test 𝒞(6) == [4,1]
    @test 𝒞(7) == [3]
    @test 𝒞(8) == [3,4]
    @test 𝒞(9) == [2]
    𝒜 = Adjacency{0}(t)
    @test 𝒜(1) == [2,5,3]
    @test 𝒜(2) == [4,6,1]
    @test 𝒜(3) == [1,5,4]
    @test 𝒜(4) == [3,6,2]
    @test 𝒜(5) == [6,3,1]
    @test 𝒜(6) == [2,4,5]

    # invalid relations
    elems = connect.([(1,2,3),(4,3,2)])
    t = HalfEdgeTopology(elems)
    @test_throws AssertionError Boundary{3,0}(t)
    @test_throws AssertionError Coboundary{0,3}(t)
    @test_throws AssertionError Adjacency{3}(t)
    @test_throws AssertionError Boundary{0,2}(t)
    @test_throws AssertionError Coboundary{2,0}(t)
  end

  @testset "SimpleTopology" begin
    elems = connect.([(1,2,3),(4,3,2)])
    t = SimpleTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1,2,3]
    @test ∂(2) == [4,3,2]
  end
end
