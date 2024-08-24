@testitem "Pointification" setup = [Setup] begin
  p = cart(0, 0)
  @test pointify(p) == [cart(0, 0)]

  sphere = Sphere(cart(0, 0), T(1))
  points = pointify(sphere)
  @test all(∈(sphere), points)

  ball = Ball(cart(0, 0), T(1))
  points = pointify(ball)
  @test all(∈(boundary(ball)), points)

  verts = [cart(0, 0), cart(1, 1)]
  segment = Segment(verts...)
  points = pointify(segment)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1)]
  triangle = Triangle(verts...)
  points = pointify(triangle)
  @test points == verts

  verts = [cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1)]
  quadrangle = Quadrangle(verts...)
  points = pointify(quadrangle)
  @test points == verts

  tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  multi = Multi([tri, quad])
  points = pointify(multi)
  @test points == [pointify(tri); pointify(quad)]

  tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  gset = GeometrySet([tri, quad])
  points = pointify(gset)
  @test points == [pointify(tri); pointify(quad)]

  pts = randpoint2(100)
  pset = PointSet(pts)
  @test pointify(pset) == pts

  grid = cartgrid(10, 10)
  @test pointify(grid) == vertices(grid)

  grid = cartgrid(10, 10)
  mesh = convert(SimpleMesh, grid)
  points = pointify(mesh)
  @test points == vertices(mesh)
  @test points == vertices(grid)
end
