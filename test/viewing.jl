@testset "Viewing" begin
  g = CartesianGrid{T}(10, 10)
  v = view(g, 1:3)
  @test unview(v) == (g, 1:3)
  @test unview(g) == (g, 1:100)

  g = CartesianGrid{T}(10, 10)
  b = Box(P2(1, 1), P2(5, 5))
  v = view(g, b)
  @test v == CartesianGrid(P2(0, 0), P2(6, 6), dims=(6, 6))

  p = PointSet(collect(vertices(g)))
  v = view(p, b)
  @test centroid(v, 1) == P2(1, 1)
  @test centroid(v, nelements(v)) == P2(5, 5)

  g = CartesianGrid{T}(10, 10)
  p = PointSet(collect(vertices(g)))
  b = Ball(P2(0, 0), T(2))
  v = view(g, b)
  @test nelements(v) == 4
  @test v[1] == g[1]
  v = view(p, b)
  @test nelements(v) == 6
  @test coordinates.(v) == V2[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)]

  # convex polygons
  tri = Triangle(P2(5, 7), P2(10, 12), P2(15, 7))
  pent = Pentagon(P2(6, 1), P2(2, 10), P2(10, 16), P2(18, 10), P2(14, 1))

  grid = CartesianGrid{T}(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[10, 10] ∈ indices(grid, tri)
  @test linds[10, 6] ∈ indices(grid, pent)

  grid = CartesianGrid(P2(-2, -2), P2(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[21, 7] ∈ indices(grid, tri)
  @test linds[21, 4] ∈ indices(grid, pent)

  grid = CartesianGrid(P2(-100, -100), P2(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, tri)
  @test linds[55, 53] ∈ indices(grid, pent)

  # non-convex polygons
  poly1 = PolyArea(P2[(3, 3), (9, 9), (3, 15), (17, 15), (17, 3)])
  poly2 = PolyArea([pointify(pent), pointify(tri)])

  grid = CartesianGrid{T}(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[12, 6] ∈ indices(grid, poly1)
  @test linds[10, 3] ∈ indices(grid, poly2)

  grid = CartesianGrid(P2(-2, -2), P2(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[22, 6] ∈ indices(grid, poly1)
  @test linds[17, 4] ∈ indices(grid, poly2)

  grid = CartesianGrid(P2(-100, -100), P2(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, poly1)
  @test linds[55, 53] ∈ indices(grid, poly2)

  # rotate
  poly1 = poly1 |> Rotate(Angle2d(π / 2))
  poly2 = poly2 |> Rotate(Angle2d(π / 2))

  grid = CartesianGrid(P2(-20, 0), P2(0, 20), T.((1, 1)))
  linds = LinearIndices(size(grid))
  @test linds[12, 12] ∈ indices(grid, poly1)
  @test linds[16, 11] ∈ indices(grid, poly2)

  grid = CartesianGrid(P2(-22, -2), P2(0, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[26, 8] ∈ indices(grid, poly1)
  @test linds[36, 9] ∈ indices(grid, poly2)

  grid = CartesianGrid(P2(-100, -100), P2(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[46, 57] ∈ indices(grid, poly1)
  @test linds[48, 55] ∈ indices(grid, poly2)

  # multi
  multi = Multi([tri, pent])
  grid = CartesianGrid{T}(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[10, 10] ∈ indices(grid, multi)
  @test linds[10, 6] ∈ indices(grid, multi)

  # clipping
  tri = Triangle(P2(-4, 10), P2(5, 19), P2(5, 1))
  grid = CartesianGrid{T}(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[3, 10] ∈ indices(grid, tri)

  # out of grid
  tri = Triangle(P2(-12, 8), P2(-8, 14), P2(-4, 8))
  grid = CartesianGrid{T}(20, 20)
  @test isempty(indices(grid, tri))

  # chain
  seg = Segment(P2(2, 12), P2(16, 18))
  rope = Rope(P2(8, 1), P2(5, 9), P2(9, 13), P2(17, 10))
  ring = Ring(P2(8, 1), P2(5, 9), P2(9, 13), P2(17, 10))
  grid = CartesianGrid{T}(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[9, 15] ∈ indices(grid, seg)
  @test linds[7, 11] ∈ indices(grid, rope)
  @test linds[12, 5] ∈ indices(grid, ring)
end
