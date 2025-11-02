@testitem "Indices" setup = [Setup] begin
  g = cartgrid(10, 10)
  b = Box(cart(1, 1), cart(5, 5))
  v = view(g, b)
  @test v == CartesianGrid(cart(0, 0), cart(6, 6), dims=(6, 6))

  p = PointSet(vertices(g))
  v = view(p, b)
  @test centroid(v, 1) == cart(1, 1)
  @test centroid(v, nelements(v)) == cart(5, 5)

  # boxes
  g = cartgrid(10, 10)
  p = PointSet(vertices(g))
  b = Ball(cart(0, 0), T(2))
  v = view(g, b)
  @test nelements(v) == 4
  @test v[1] == g[1]
  v = view(p, b)
  @test nelements(v) == 6
  @test to.(v) == vector.([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)])

  b1 = Box(cart(0.6, 0.7), cart(1.0, 1.0))
  b2 = Box(cart(0.6, 0.7), cart(1.2, 1.2))
  b3 = Box(cart(0.0, 0.0), cart(0.6, 0.3))
  b4 = Box(cart(-0.2, -0.2), cart(0.6, 0.3))
  b5 = Box(cart(0.25, 0.15), cart(0.95, 0.85))
  b6 = Box(cart(0.35, 0.25), cart(0.85, 0.75))
  b7 = Box(cart(0.05, 0.05), cart(0.65, 0.35))
  b8 = Box(cart(0.55, 0.65), cart(0.95, 0.95))
  x = range(zero(T), stop=one(T), length=6)
  y = T[0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
  g = RectilinearGrid(x, y)
  linds = LinearIndices(size(g))
  @test issetequal(indices(g, b1), [linds[4, 4], linds[5, 4], linds[4, 5], linds[5, 5]])
  @test issetequal(indices(g, b2), [linds[4, 4], linds[5, 4], linds[4, 5], linds[5, 5]])
  @test issetequal(indices(g, b3), [linds[1, 1], linds[2, 1], linds[3, 1], linds[1, 2], linds[2, 2], linds[3, 2]])
  @test issetequal(indices(g, b4), [linds[1, 1], linds[2, 1], linds[3, 1], linds[1, 2], linds[2, 2], linds[3, 2]])
  @test issetequal(
    indices(g, b5),
    [
      linds[2, 2],
      linds[3, 2],
      linds[4, 2],
      linds[5, 2],
      linds[2, 3],
      linds[3, 3],
      linds[4, 3],
      linds[5, 3],
      linds[2, 4],
      linds[3, 4],
      linds[4, 4],
      linds[5, 4]
    ]
  )
  @test issetequal(
    indices(g, b6),
    [
      linds[2, 2],
      linds[3, 2],
      linds[4, 2],
      linds[5, 2],
      linds[2, 3],
      linds[3, 3],
      linds[4, 3],
      linds[5, 3],
      linds[2, 4],
      linds[3, 4],
      linds[4, 4],
      linds[5, 4]
    ]
  )
  @test issetequal(
    indices(g, b7),
    [
      linds[1, 1],
      linds[2, 1],
      linds[3, 1],
      linds[4, 1],
      linds[1, 2],
      linds[2, 2],
      linds[3, 2],
      linds[4, 2],
      linds[1, 3],
      linds[2, 3],
      linds[3, 3],
      linds[4, 3]
    ]
  )
  @test issetequal(
    indices(g, b8),
    [
      linds[3, 3],
      linds[4, 3],
      linds[5, 3],
      linds[3, 4],
      linds[4, 4],
      linds[5, 4],
      linds[3, 5],
      linds[4, 5],
      linds[5, 5]
    ]
  )

  # convex polygons
  tri = Triangle(cart(5, 7), cart(10, 12), cart(15, 7))
  pent = Pentagon(cart(6, 1), cart(2, 10), cart(10, 16), cart(18, 10), cart(14, 1))

  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[10, 10] ∈ indices(grid, tri)
  @test linds[10, 6] ∈ indices(grid, pent)

  grid = CartesianGrid(cart(-2, -2), cart(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[21, 7] ∈ indices(grid, tri)
  @test linds[21, 4] ∈ indices(grid, pent)

  grid = CartesianGrid(cart(-100, -100), cart(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, tri)
  @test linds[55, 53] ∈ indices(grid, pent)

  # non-convex polygons
  poly1 = PolyArea(cart.([(3, 3), (9, 9), (3, 15), (17, 15), (17, 3)]))
  poly2 = PolyArea([pointify(pent), pointify(tri)])

  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[12, 6] ∈ indices(grid, poly1)
  @test linds[10, 3] ∈ indices(grid, poly2)

  grid = CartesianGrid(cart(-2, -2), cart(20, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[22, 6] ∈ indices(grid, poly1)
  @test linds[17, 4] ∈ indices(grid, poly2)

  grid = CartesianGrid(cart(-100, -100), cart(20, 20), T.((2, 2)))
  linds = LinearIndices(size(grid))
  @test linds[57, 54] ∈ indices(grid, poly1)
  @test linds[55, 53] ∈ indices(grid, poly2)

  # rotate
  poly1 = poly1 |> Rotate(Angle2d(T(π / 2)))
  poly2 = poly2 |> Rotate(Angle2d(T(π / 2)))

  grid = CartesianGrid(cart(-20, 0), cart(0, 20), T.((1, 1)))
  linds = LinearIndices(size(grid))
  @test linds[12, 12] ∈ indices(grid, poly1)
  @test linds[16, 11] ∈ indices(grid, poly2)

  grid = CartesianGrid(cart(-22, -2), cart(0, 20), T.((0.5, 1.5)))
  linds = LinearIndices(size(grid))
  @test linds[26, 8] ∈ indices(grid, poly1)
  @test linds[36, 9] ∈ indices(grid, poly2)

  grid = CartesianGrid(cart(-100, -100), cart(20, 20), T.((2, 2)))
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
  tri = Triangle(cart(-4, 10), cart(5, 19), cart(5, 1))
  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[3, 10] ∈ indices(grid, tri)

  # out of grid
  tri = Triangle(cart(-12, 8), cart(-8, 14), cart(-4, 8))
  grid = cartgrid(20, 20)
  @test isempty(indices(grid, tri))

  # chain
  seg = Segment(cart(2, 12), cart(16, 18))
  rope = Rope(cart(8, 1), cart(5, 9), cart(9, 13), cart(17, 10))
  ring = Ring(cart(8, 1), cart(5, 9), cart(9, 13), cart(17, 10))
  grid = cartgrid(20, 20)
  linds = LinearIndices(size(grid))
  @test linds[9, 15] ∈ indices(grid, seg)
  @test linds[7, 11] ∈ indices(grid, rope)
  @test linds[12, 5] ∈ indices(grid, ring)

  # points
  p1 = cart(0, 0)
  p2 = cart(0.5, 0.5)
  p3 = cart(1, 1)
  p4 = cart(2, 2)
  p5 = cart(10, 10)
  p6 = cart(11, 11)
  grid = cartgrid(10, 10)
  linds = LinearIndices(size(grid))
  @test linds[1, 1] == only(indices(grid, p1))
  @test linds[1, 1] == only(indices(grid, p2))
  @test linds[1, 1] == only(indices(grid, p3))
  @test linds[2, 2] == only(indices(grid, p4))
  @test linds[10, 10] == only(indices(grid, p5))
  @test isempty(indices(grid, p6))

  # transformed grid
  g = cartgrid(20, 20)
  poly1 = PolyArea(cart.([(3, 3), (9, 9), (3, 15), (17, 15), (17, 3)]))
  expectedinds = indices(g, poly1)
  linds = LinearIndices(size(g))
  @test linds[12, 6] ∈ expectedinds
  taff = Rotate(Angle2d(T(π / 2)))
  tmorph = Morphological(c -> Cartesian(c.x + c.y, c.y, c.x - c.y))
  tseq = taff → tmorph
  grid = TransformedGrid(g, tseq)
  poly2 = taff(poly1)
  actualinds = indices(grid, poly2)
  @test actualinds == expectedinds
end
