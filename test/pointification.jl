@testset "Pointification" begin
  point = P2(0, 0)
  @test pointify(point) == [P2(0, 0)]

  sphere = Sphere(P2(0, 0), T(1))
  points = pointify(sphere)
  @test all(∈(sphere), points)

  ball = Ball(P2(0, 0), T(1))
  points = pointify(ball)
  @test all(∈(boundary(ball)), points)

  verts = [P2(0, 0), P2(1, 1)]
  segment = Segment(verts...)
  points = pointify(segment)
  @test points == verts

  verts = [P2(0, 0), P2(1, 0), P2(1, 1)]
  triangle = Triangle(verts...)
  points = pointify(triangle)
  @test points == verts

  verts = [P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1)]
  quadrangle = Quadrangle(verts...)
  points = pointify(quadrangle)
  @test points == verts

  tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
  quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
  multi = Multi([tri, quad])
  points = pointify(multi)
  @test points == [pointify(tri); pointify(quad)]

  tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
  quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
  gset = GeometrySet([tri, quad])
  points = pointify(gset)
  @test points == [pointify(tri); pointify(quad)]

  pts = rand(P2, 100)
  pset = PointSet(pts)
  @test pointify(pset) == pts

  grid = CartesianGrid{T}(10, 10)
  @test pointify(grid) == vertices(grid)

  grid = CartesianGrid{T}(10, 10)
  mesh = convert(SimpleMesh, grid)
  points = pointify(mesh)
  @test points == vertices(mesh)
  @test points == vertices(grid)
end
