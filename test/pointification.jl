@testset "Pointification" begin
  p = point(0, 0)
  @test pointify(p) == [point(0, 0)]

  sphere = Sphere(point(0, 0), T(1))
  points = pointify(sphere)
  @test all(∈(sphere), points)

  ball = Ball(point(0, 0), T(1))
  points = pointify(ball)
  @test all(∈(boundary(ball)), points)

  verts = [point(0, 0), point(1, 1)]
  segment = Segment(verts...)
  points = pointify(segment)
  @test points == verts

  verts = [point(0, 0), point(1, 0), point(1, 1)]
  triangle = Triangle(verts...)
  points = pointify(triangle)
  @test points == verts

  verts = [point(0, 0), point(1, 0), point(1, 1), point(0, 1)]
  quadrangle = Quadrangle(verts...)
  points = pointify(quadrangle)
  @test points == verts

  tri = Triangle(point(0, 0), point(1, 0), point(1, 1))
  quad = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
  multi = Multi([tri, quad])
  points = pointify(multi)
  @test points == [pointify(tri); pointify(quad)]

  tri = Triangle(point(0, 0), point(1, 0), point(1, 1))
  quad = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
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
