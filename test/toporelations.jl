@testset "TopologicalRelation" begin
  @testset "GridTopology" begin
    # 3 segments
    t = GridTopology(3)
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1, 2]
    @test ∂(2) == [2, 3]
    @test ∂(3) == [3, 4]

    # quadrangles in 2D grid
    t = GridTopology(2, 3)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1, 2, 5, 4]
    @test ∂(2) == [2, 3, 6, 5]
    @test ∂(3) == [4, 5, 8, 7]
    @test ∂(4) == [5, 6, 9, 8]
    @test ∂(5) == [7, 8, 11, 10]
    @test ∂(6) == [8, 9, 12, 11]

    # segments of quadrangles in 2D grid
    t = GridTopology(2, 3)
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1, 4]
    @test ∂(2) == [2, 5]
    @test ∂(3) == [3, 6]
    @test ∂(4) == [4, 7]
    @test ∂(5) == [5, 8]
    @test ∂(6) == [6, 9]
    @test ∂(7) == [7, 10]
    @test ∂(8) == [8, 11]
    @test ∂(9) == [9, 12]
    @test ∂(10) == [1, 2]
    @test ∂(11) == [4, 5]
    @test ∂(12) == [7, 8]
    @test ∂(13) == [10, 11]
    @test ∂(14) == [2, 3]
    @test ∂(15) == [5, 6]
    @test ∂(16) == [8, 9]
    @test ∂(17) == [11, 12]

    # segments of quadrangles in 2D (periodic) grid
    t = GridTopology((2, 2), (true, false))
    ∂ = Boundary{1,0}(t)
    @test nfacets(t) == 10
    @test ∂(1) == [1, 3]
    @test ∂(2) == [2, 4]
    @test ∂(3) == [3, 5]
    @test ∂(4) == [4, 6]
    @test ∂(5) == [1, 2]
    @test ∂(6) == [3, 4]
    @test ∂(7) == [5, 6]
    @test ∂(8) == [2, 1]
    @test ∂(9) == [4, 3]
    @test ∂(10) == [6, 5]

    # segments of quadrangles in 2D (periodic) grid
    t = GridTopology((2, 2), (false, true))
    ∂ = Boundary{1,0}(t)
    @test nfacets(t) == 10
    @test ∂(1) == [1, 4]
    @test ∂(2) == [2, 5]
    @test ∂(3) == [3, 6]
    @test ∂(4) == [4, 1]
    @test ∂(5) == [5, 2]
    @test ∂(6) == [6, 3]
    @test ∂(7) == [1, 2]
    @test ∂(8) == [4, 5]
    @test ∂(9) == [2, 3]
    @test ∂(10) == [5, 6]

    # segments of quadrangles in 2D (periodic) grid
    t = GridTopology((2, 2), (true, true))
    ∂ = Boundary{1,0}(t)
    @test nfacets(t) == 8
    @test ∂(1) == [1, 3]
    @test ∂(2) == [2, 4]
    @test ∂(3) == [3, 1]
    @test ∂(4) == [4, 2]
    @test ∂(5) == [1, 2]
    @test ∂(6) == [3, 4]
    @test ∂(7) == [2, 1]
    @test ∂(8) == [4, 3]

    # quadrangles of hexahedrons in 3D grid
    t = GridTopology(2, 2, 2)
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 13, 14, 25, 26]
    @test ∂(2) == [2, 3, 16, 17, 28, 29]
    @test ∂(3) == [4, 5, 14, 15, 31, 32]
    @test ∂(4) == [5, 6, 17, 18, 34, 35]
    @test ∂(5) == [7, 8, 19, 20, 26, 27]
    @test ∂(6) == [8, 9, 22, 23, 29, 30]
    @test ∂(7) == [10, 11, 20, 21, 32, 33]
    @test ∂(8) == [11, 12, 23, 24, 35, 36]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (true, false, false))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 9, 10, 21, 22]
    @test ∂(2) == [2, 1, 12, 13, 24, 25]
    @test ∂(3) == [3, 4, 10, 11, 27, 28]
    @test ∂(4) == [4, 3, 13, 14, 30, 31]
    @test ∂(5) == [5, 6, 15, 16, 22, 23]
    @test ∂(6) == [6, 5, 18, 19, 25, 26]
    @test ∂(7) == [7, 8, 16, 17, 28, 29]
    @test ∂(8) == [8, 7, 19, 20, 31, 32]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (false, true, false))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 13, 14, 21, 22]
    @test ∂(2) == [2, 3, 15, 16, 24, 25]
    @test ∂(3) == [4, 5, 14, 13, 27, 28]
    @test ∂(4) == [5, 6, 16, 15, 30, 31]
    @test ∂(5) == [7, 8, 17, 18, 22, 23]
    @test ∂(6) == [8, 9, 19, 20, 25, 26]
    @test ∂(7) == [10, 11, 18, 17, 28, 29]
    @test ∂(8) == [11, 12, 20, 19, 31, 32]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (false, false, true))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 13, 14, 25, 26]
    @test ∂(2) == [2, 3, 16, 17, 27, 28]
    @test ∂(3) == [4, 5, 14, 15, 29, 30]
    @test ∂(4) == [5, 6, 17, 18, 31, 32]
    @test ∂(5) == [7, 8, 19, 20, 26, 25]
    @test ∂(6) == [8, 9, 22, 23, 28, 27]
    @test ∂(7) == [10, 11, 20, 21, 30, 29]
    @test ∂(8) == [11, 12, 23, 24, 32, 31]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (true, true, false))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 9, 10, 17, 18]
    @test ∂(2) == [2, 1, 11, 12, 20, 21]
    @test ∂(3) == [3, 4, 10, 9, 23, 24]
    @test ∂(4) == [4, 3, 12, 11, 26, 27]
    @test ∂(5) == [5, 6, 13, 14, 18, 19]
    @test ∂(6) == [6, 5, 15, 16, 21, 22]
    @test ∂(7) == [7, 8, 14, 13, 24, 25]
    @test ∂(8) == [8, 7, 16, 15, 27, 28]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (true, false, true))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 9, 10, 21, 22]
    @test ∂(2) == [2, 1, 12, 13, 23, 24]
    @test ∂(3) == [3, 4, 10, 11, 25, 26]
    @test ∂(4) == [4, 3, 13, 14, 27, 28]
    @test ∂(5) == [5, 6, 15, 16, 22, 21]
    @test ∂(6) == [6, 5, 18, 19, 24, 23]
    @test ∂(7) == [7, 8, 16, 17, 26, 25]
    @test ∂(8) == [8, 7, 19, 20, 28, 27]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (false, true, true))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 13, 14, 21, 22]
    @test ∂(2) == [2, 3, 15, 16, 23, 24]
    @test ∂(3) == [4, 5, 14, 13, 25, 26]
    @test ∂(4) == [5, 6, 16, 15, 27, 28]
    @test ∂(5) == [7, 8, 17, 18, 22, 21]
    @test ∂(6) == [8, 9, 19, 20, 24, 23]
    @test ∂(7) == [10, 11, 18, 17, 26, 25]
    @test ∂(8) == [11, 12, 20, 19, 28, 27]

    # quadrangles of hexahedrons in 3D (periodic) grid
    t = GridTopology((2, 2, 2), (true, true, true))
    ∂ = Boundary{3,2}(t)
    @test ∂(1) == [1, 2, 9, 10, 17, 18]
    @test ∂(2) == [2, 1, 11, 12, 19, 20]
    @test ∂(3) == [3, 4, 10, 9, 21, 22]
    @test ∂(4) == [4, 3, 12, 11, 23, 24]
    @test ∂(5) == [5, 6, 13, 14, 18, 17]
    @test ∂(6) == [6, 5, 15, 16, 20, 19]
    @test ∂(7) == [7, 8, 14, 13, 22, 21]
    @test ∂(8) == [8, 7, 16, 15, 24, 23]

    # edges of quadrangles in 2D grid
    t = GridTopology(3, 4)
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 2, 17, 18]
    @test ∂(2) == [2, 3, 22, 23]
    @test ∂(3) == [3, 4, 27, 28]
    @test ∂(4) == [5, 6, 18, 19]
    @test ∂(5) == [6, 7, 23, 24]
    @test ∂(6) == [7, 8, 28, 29]
    @test ∂(7) == [9, 10, 19, 20]
    @test ∂(8) == [10, 11, 24, 25]
    @test ∂(9) == [11, 12, 29, 30]
    @test ∂(10) == [13, 14, 20, 21]
    @test ∂(11) == [14, 15, 25, 26]
    @test ∂(12) == [15, 16, 30, 31]

    # edges of quadrangles in 2D (periodic) grid
    t = GridTopology((3, 4), (true, false))
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 2, 13, 14]
    @test ∂(2) == [2, 3, 18, 19]
    @test ∂(3) == [3, 1, 23, 24]
    @test ∂(4) == [4, 5, 14, 15]
    @test ∂(5) == [5, 6, 19, 20]
    @test ∂(6) == [6, 4, 24, 25]
    @test ∂(7) == [7, 8, 15, 16]
    @test ∂(8) == [8, 9, 20, 21]
    @test ∂(9) == [9, 7, 25, 26]
    @test ∂(10) == [10, 11, 16, 17]
    @test ∂(11) == [11, 12, 21, 22]
    @test ∂(12) == [12, 10, 26, 27]

    # edges of quadrangles in 2D (periodic) grid
    t = GridTopology((3, 4), (false, true))
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 2, 17, 18]
    @test ∂(2) == [2, 3, 21, 22]
    @test ∂(3) == [3, 4, 25, 26]
    @test ∂(4) == [5, 6, 18, 19]
    @test ∂(5) == [6, 7, 22, 23]
    @test ∂(6) == [7, 8, 26, 27]
    @test ∂(7) == [9, 10, 19, 20]
    @test ∂(8) == [10, 11, 23, 24]
    @test ∂(9) == [11, 12, 27, 28]
    @test ∂(10) == [13, 14, 20, 17]
    @test ∂(11) == [14, 15, 24, 21]
    @test ∂(12) == [15, 16, 28, 25]

    # edges of quadrangles in 2D (periodic) grid
    t = GridTopology((3, 4), (true, true))
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 2, 13, 14]
    @test ∂(2) == [2, 3, 17, 18]
    @test ∂(3) == [3, 1, 21, 22]
    @test ∂(4) == [4, 5, 14, 15]
    @test ∂(5) == [5, 6, 18, 19]
    @test ∂(6) == [6, 4, 22, 23]
    @test ∂(7) == [7, 8, 15, 16]
    @test ∂(8) == [8, 9, 19, 20]
    @test ∂(9) == [9, 7, 23, 24]
    @test ∂(10) == [10, 11, 16, 13]
    @test ∂(11) == [11, 12, 20, 17]
    @test ∂(12) == [12, 10, 24, 21]

    # 2x3x2 hexahedrons
    t = GridTopology(2, 3, 2)
    ∂ = Boundary{3,0}(t)
    @test ∂(1) == [1, 2, 5, 4, 13, 14, 17, 16]
    @test ∂(2) == [2, 3, 6, 5, 14, 15, 18, 17]
    @test ∂(3) == [4, 5, 8, 7, 16, 17, 20, 19]
    @test ∂(12) == [20, 21, 24, 23, 32, 33, 36, 35]

    # quadrangles in 2D grid
    t = GridTopology(2, 3)
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [2, 3]
    @test 𝒜(2) == [1, 4]
    @test 𝒜(3) == [4, 1, 5]
    @test 𝒜(4) == [3, 2, 6]
    @test 𝒜(5) == [6, 3]
    @test 𝒜(6) == [5, 4]

    # quadrangles in 2D grid
    t = GridTopology(3, 3)
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [2, 4]
    @test 𝒜(2) == [1, 3, 5]
    @test 𝒜(3) == [2, 6]
    @test 𝒜(4) == [5, 1, 7]
    @test 𝒜(5) == [4, 6, 2, 8]
    @test 𝒜(6) == [5, 3, 9]
    @test 𝒜(7) == [8, 4]
    @test 𝒜(8) == [7, 9, 5]
    @test 𝒜(9) == [8, 6]

    # quadrangles in 2D grid with periodicity
    t = GridTopology((3, 3), (true, false))
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [3, 2, 4]
    @test 𝒜(2) == [1, 3, 5]
    @test 𝒜(3) == [2, 1, 6]
    @test 𝒜(4) == [6, 5, 1, 7]
    @test 𝒜(5) == [4, 6, 2, 8]
    @test 𝒜(6) == [5, 4, 3, 9]
    @test 𝒜(7) == [9, 8, 4]
    @test 𝒜(8) == [7, 9, 5]
    @test 𝒜(9) == [8, 7, 6]

    # quadrangles in 2D grid with periodicity
    t = GridTopology((3, 3), (true, true))
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [3, 2, 7, 4]
    @test 𝒜(2) == [1, 3, 8, 5]
    @test 𝒜(3) == [2, 1, 9, 6]
    @test 𝒜(4) == [6, 5, 1, 7]
    @test 𝒜(5) == [4, 6, 2, 8]
    @test 𝒜(6) == [5, 4, 3, 9]
    @test 𝒜(7) == [9, 8, 4, 1]
    @test 𝒜(8) == [7, 9, 5, 2]
    @test 𝒜(9) == [8, 7, 6, 3]

    # quadrangles in 3D grid
    t = GridTopology(2, 2, 2)
    𝒜 = Adjacency{3}(t)
    @test 𝒜(1) == [2, 3, 5]
    @test 𝒜(2) == [1, 4, 6]
    @test 𝒜(3) == [4, 1, 7]
    @test 𝒜(4) == [3, 2, 8]
    @test 𝒜(5) == [6, 7, 1]
    @test 𝒜(6) == [5, 8, 2]
    @test 𝒜(7) == [8, 5, 3]
    @test 𝒜(8) == [7, 6, 4]

    # quadrangles in 3D grid
    t = GridTopology(3, 2, 2)
    𝒜 = Adjacency{3}(t)
    @test 𝒜(1) == [2, 4, 7]
    @test 𝒜(2) == [1, 3, 5, 8]
    @test 𝒜(3) == [2, 6, 9]
    @test 𝒜(4) == [5, 1, 10]
    @test 𝒜(5) == [4, 6, 2, 11]
    @test 𝒜(6) == [5, 3, 12]
    @test 𝒜(7) == [8, 10, 1]
    @test 𝒜(8) == [7, 9, 11, 2]
    @test 𝒜(9) == [8, 12, 3]
    @test 𝒜(10) == [11, 7, 4]
    @test 𝒜(11) == [10, 12, 8, 5]
    @test 𝒜(12) == [11, 9, 6]

    # quadrangles in 3D grid with periodicity
    t = GridTopology((3, 2, 2), (true, false, false))
    𝒜 = Adjacency{3}(t)
    @test 𝒜(1) == [3, 2, 4, 7]
    @test 𝒜(2) == [1, 3, 5, 8]
    @test 𝒜(3) == [2, 1, 6, 9]
    @test 𝒜(4) == [6, 5, 1, 10]
    @test 𝒜(5) == [4, 6, 2, 11]
    @test 𝒜(6) == [5, 4, 3, 12]
    @test 𝒜(7) == [9, 8, 10, 1]
    @test 𝒜(8) == [7, 9, 11, 2]
    @test 𝒜(9) == [8, 7, 12, 3]
    @test 𝒜(10) == [12, 11, 7, 4]
    @test 𝒜(11) == [10, 12, 8, 5]
    @test 𝒜(12) == [11, 10, 9, 6]

    # vertices in 2D grid
    t = GridTopology(2, 2)
    𝒜 = Adjacency{0}(t)
    @test 𝒜(1) == [2, 4]
    @test 𝒜(2) == [1, 3, 5]
    @test 𝒜(3) == [2, 6]
    @test 𝒜(4) == [5, 1, 7]
    @test 𝒜(5) == [4, 6, 2, 8]
    @test 𝒜(6) == [5, 3, 9]
    @test 𝒜(7) == [8, 4]
    @test 𝒜(8) == [7, 9, 5]
    @test 𝒜(9) == [8, 6]

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
    elems = connect.([(1, 2, 3), (4, 3, 2)])
    t = HalfEdgeTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [2, 3, 1]
    @test ∂(2) == [3, 2, 4]
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 3, 2]
    @test ∂(2) == [1, 4, 5]
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [3, 2]
    @test ∂(2) == [1, 2]
    @test ∂(3) == [3, 1]
    @test ∂(4) == [2, 4]
    @test ∂(5) == [4, 3]
    𝒞 = Coboundary{0,1}(t)
    @test 𝒞(1) == [2, 3]
    @test 𝒞(2) == [4, 1, 2]
    @test 𝒞(3) == [3, 1, 5]
    @test 𝒞(4) == [5, 4]
    𝒞 = Coboundary{0,2}(t)
    @test 𝒞(1) == [1]
    @test 𝒞(2) == [2, 1]
    @test 𝒞(3) == [1, 2]
    @test 𝒞(4) == [2]
    𝒞 = Coboundary{1,2}(t)
    @test 𝒞(1) == [2, 1]
    @test 𝒞(2) == [1]
    @test 𝒞(3) == [1]
    @test 𝒞(4) == [2]
    @test 𝒞(5) == [2]
    𝒜 = Adjacency{0}(t)
    @test 𝒜(1) == [2, 3]
    @test 𝒜(2) == [4, 3, 1]
    @test 𝒜(3) == [1, 2, 4]
    @test 𝒜(4) == [3, 2]

    # 2 triangles + 2 quadrangles
    elems = connect.([(1, 2, 6, 5), (2, 4, 6), (4, 3, 5, 6), (1, 5, 3)])
    t = HalfEdgeTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1, 2, 6, 5]
    @test ∂(2) == [6, 2, 4]
    @test ∂(3) == [6, 4, 3, 5]
    @test ∂(4) == [3, 1, 5]
    ∂ = Boundary{2,1}(t)
    @test ∂(1) == [1, 3, 5, 6]
    @test ∂(2) == [3, 9, 4]
    @test ∂(3) == [4, 7, 8, 5]
    @test ∂(4) == [2, 6, 8]
    ∂ = Boundary{1,0}(t)
    @test ∂(1) == [1, 2]
    @test ∂(2) == [3, 1]
    @test ∂(3) == [6, 2]
    @test ∂(4) == [4, 6]
    @test ∂(5) == [5, 6]
    @test ∂(6) == [1, 5]
    @test ∂(7) == [4, 3]
    @test ∂(8) == [3, 5]
    @test ∂(9) == [2, 4]
    𝒞 = Coboundary{0,1}(t)
    @test 𝒞(1) == [1, 6, 2]
    @test 𝒞(2) == [9, 3, 1]
    @test 𝒞(3) == [2, 8, 7]
    @test 𝒞(4) == [7, 4, 9]
    @test 𝒞(5) == [5, 8, 6]
    @test 𝒞(6) == [3, 4, 5]
    𝒞 = Coboundary{0,2}(t)
    @test 𝒞(1) == [1, 4]
    @test 𝒞(2) == [2, 1]
    @test 𝒞(3) == [4, 3]
    @test 𝒞(4) == [3, 2]
    @test 𝒞(5) == [3, 4, 1]
    @test 𝒞(6) == [2, 3, 1]
    𝒞 = Coboundary{1,2}(t)
    @test 𝒞(1) == [1]
    @test 𝒞(2) == [4]
    @test 𝒞(3) == [2, 1]
    @test 𝒞(4) == [2, 3]
    @test 𝒞(5) == [3, 1]
    @test 𝒞(6) == [4, 1]
    @test 𝒞(7) == [3]
    @test 𝒞(8) == [3, 4]
    @test 𝒞(9) == [2]
    𝒜 = Adjacency{0}(t)
    @test 𝒜(1) == [2, 5, 3]
    @test 𝒜(2) == [4, 6, 1]
    @test 𝒜(3) == [1, 5, 4]
    @test 𝒜(4) == [3, 6, 2]
    @test 𝒜(5) == [6, 3, 1]
    @test 𝒜(6) == [2, 4, 5]

    # 2 triangles + 2 quadrangles
    elems = connect.([(1, 2, 6, 5), (2, 4, 6), (4, 3, 5, 6), (1, 5, 3)])
    t = HalfEdgeTopology(elems)
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [2, 3, 4]
    @test 𝒜(2) == [1, 3]
    @test 𝒜(3) == [2, 4, 1]
    @test 𝒜(4) == [1, 3]

    # 4 quadrangles in a grid
    elems = connect.([(1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8, 7), (5, 6, 9, 8)])
    t = HalfEdgeTopology(elems)
    𝒜 = Adjacency{2}(t)
    @test 𝒜(1) == [3, 2]
    @test 𝒜(2) == [1, 4]
    @test 𝒜(3) == [1, 4]
    @test 𝒜(4) == [3, 2]

    # invalid relations
    elems = connect.([(1, 2, 3), (4, 3, 2)])
    t = HalfEdgeTopology(elems)
    @test_throws AssertionError Boundary{3,0}(t)
    @test_throws AssertionError Coboundary{0,3}(t)
    @test_throws AssertionError Adjacency{3}(t)
    @test_throws AssertionError Boundary{0,2}(t)
    @test_throws AssertionError Coboundary{2,0}(t)
  end

  @testset "SimpleTopology" begin
    elems = connect.([(1, 2, 3), (4, 3, 2)])
    t = SimpleTopology(elems)
    ∂ = Boundary{2,0}(t)
    @test ∂(1) == [1, 2, 3]
    @test ∂(2) == [4, 3, 2]
  end
end
