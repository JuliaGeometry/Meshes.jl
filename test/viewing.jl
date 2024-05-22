@testset "Viewing" begin
  g = cartgrid(10, 10)
  v = view(g, 1:3)
  @test parent(v) == g
  @test parentindices(v) == 1:3
  @test parent(g) == g
  @test parentindices(g) == 1:100

  g = cartgrid(10, 10)
  b = Box(point(1, 1), point(5, 5))
  v = view(g, b)
  @test v == CartesianGrid(point(0, 0), point(6, 6), dims=(6, 6))

  p = PointSet(collect(vertices(g)))
  v = view(p, b)
  @test centroid(v, 1) == point(1, 1)
  @test centroid(v, nelements(v)) == point(5, 5)

  g = cartgrid(10, 10)
  p = PointSet(collect(vertices(g)))
  b = Ball(point(0, 0), T(2))
  v = view(g, b)
  @test nelements(v) == 4
  @test v[1] == g[1]
  v = view(p, b)
  @test nelements(v) == 6
  @test to.(v) == vector.([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)])

  # convex polygons
  tri = Triangle(point(5, 7), point(10, 12), point(15, 7))
  pent = Pentagon(point(6, 1), point(2, 10), point(10, 16), point(18, 10), point(14, 1))

  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[10, 10] ∈ indices(grid, tri)
  @test linds[10, 6] ∈ indices(grid, pent)

  grid = CartesianGrid(point(-2, -2), point(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[21, 7] ∈ indices(grid, tri)
  @test linds[21, 4] ∈ indices(grid, pent)

  grid = CartesianGrid(point(-100, -100), point(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, tri)
  @test linds[55, 53] ∈ indices(grid, pent)

  # non-convex polygons
  poly1 = PolyArea(point.([(3, 3), (9, 9), (3, 15), (17, 15), (17, 3)]))
  poly2 = PolyArea([pointify(pent), pointify(tri)])

  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[12, 6] ∈ indices(grid, poly1)
  @test linds[10, 3] ∈ indices(grid, poly2)

  grid = CartesianGrid(point(-2, -2), point(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[22, 6] ∈ indices(grid, poly1)
  @test linds[17, 4] ∈ indices(grid, poly2)

  grid = CartesianGrid(point(-100, -100), point(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, poly1)
  @test linds[55, 53] ∈ indices(grid, poly2)

  # rotate
  poly1 = poly1 |> Rotate(Angle2d(T(π / 2)))
  poly2 = poly2 |> Rotate(Angle2d(T(π / 2)))

  grid = CartesianGrid(point(-20, 0), point(0, 20), T.((1, 1)))
  linds = LinearIndices(size(grid))
  @test linds[12, 12] ∈ indices(grid, poly1)
  @test linds[16, 11] ∈ indices(grid, poly2)

  grid = CartesianGrid(point(-22, -2), point(0, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[26, 8] ∈ indices(grid, poly1)
  @test linds[36, 9] ∈ indices(grid, poly2)

  grid = CartesianGrid(point(-100, -100), point(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[46, 57] ∈ indices(grid, poly1)
  @test linds[48, 55] ∈ indices(grid, poly2)

  # multi
  multi = Multi([tri, pent])
  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[10, 10] ∈ indices(grid, multi)
  @test linds[10, 6] ∈ indices(grid, multi)

  # clipping
  tri = Triangle(point(-4, 10), point(5, 19), point(5, 1))
  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[3, 10] ∈ indices(grid, tri)

  # out of grid
  tri = Triangle(point(-12, 8), point(-8, 14), point(-4, 8))
  grid = cartgrid(20, 20)
  @test isempty(indices(grid, tri))

  # chain
  seg = Segment(point(2, 12), point(16, 18))
  rope = Rope(point(8, 1), point(5, 9), point(9, 13), point(17, 10))
  ring = Ring(point(8, 1), point(5, 9), point(9, 13), point(17, 10))
  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[9, 15] ∈ indices(grid, seg)
  @test linds[7, 11] ∈ indices(grid, rope)
  @test linds[12, 5] ∈ indices(grid, ring)

  # points
  p1 = point(0, 0)
  p2 = point(0.5, 0.5)
  p3 = point(1, 1)
  p4 = point(2, 2)
  p5 = point(10, 10)
  p6 = point(11, 11)
  grid = cartgrid(10, 10)
  linds = LinearIndices(size(grid))
  @test linds[1, 1] == only(indices(grid, p1))
  @test linds[1, 1] == only(indices(grid, p2))
  @test linds[1, 1] == only(indices(grid, p3))
  @test linds[2, 2] == only(indices(grid, p4))
  @test linds[10, 10] == only(indices(grid, p5))
  @test isempty(indices(grid, p6))
end
