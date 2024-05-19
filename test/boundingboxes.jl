@testset "Bounding boxes" begin
  p = point(0, 0)
  @test boundingbox(p) == Box(p, p)
  @test @allocated(boundingbox(p)) < 50

  b = Box(point(0, 0), point(1, 1))
  @test boundingbox(b) == b
  @test @allocated(boundingbox(b)) < 50

  r = Ray(point(0, 0), vector(1, 0))
  @test boundingbox(r) == Box(point(0, 0), point(T(Inf), 0))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(point(1, 1), vector(0, 1))
  @test boundingbox(r) == Box(point(1, 1), point(1, T(Inf)))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(point(1, 1), vector(-1, -1))
  @test boundingbox(r) == Box(point(T(-Inf), T(-Inf)), point(1, 1))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(point(-1, 1), vector(1, -1))
  @test boundingbox(r) == Box(point(-1, T(-Inf)), point(T(Inf), 1))
  @test @allocated(boundingbox(r)) < 50

  b = Ball(point(0, 0), T(1))
  @test boundingbox(b) == Box(point(-1, -1), point(1, 1))
  @test @allocated(boundingbox(b)) < 50
  b = Ball(point(1, 1), T(1))
  @test boundingbox(b) == Box(point(0, 0), point(2, 2))
  @test @allocated(boundingbox(b)) < 50

  s = Sphere(point(0, 0), T(1))
  @test boundingbox(s) == Box(point(-1, -1), point(1, 1))
  @test @allocated(boundingbox(s)) < 50
  s = Sphere(point(1, 1), T(1))
  @test boundingbox(s) == Box(point(0, 0), point(2, 2))
  @test @allocated(boundingbox(s)) < 50

  c = Cylinder(T(1))
  b = boundingbox(c)
  @test b == Box(point(-1, -1, 0), point(1, 1, 1))

  c = CylinderSurface(T(1))
  b = boundingbox(c)
  @test b == Box(point(-1, -1, 0), point(1, 1, 1))

  c = Cone(Disk(Plane(point(0, 0, 0), vector(0, 0, 1)), T(1)), point(0, 0, 1))
  b = boundingbox(c)
  @test b == Box(point(-1, -1, 0), point(1, 1, 1))

  c = ConeSurface(Disk(Plane(point(0, 0, 0), vector(0, 0, 1)), T(1)), point(0, 0, 1))
  b = boundingbox(c)
  @test b == Box(point(-1, -1, 0), point(1, 1, 1))

  b = Box(point(-3, -1), point(0.5, 0.5))
  s = Sphere(point(0, 0), T(2))
  m = Multi([b, s])
  d = GeometrySet([b, s])
  @test boundingbox(m) == Box(point(-3, -2), point(2, 2))
  @test boundingbox(d) == Box(point(-3, -2), point(2, 2))
  @test @allocated(boundingbox(m)) < 2500
  @test @allocated(boundingbox(d)) < 2500

  b1 = Box(point(0, 0), point(1, 1))
  b2 = Box(point(-1, -1), point(0.5, 0.5))
  m = Multi([b1, b2])
  d = GeometrySet([b1, b2])
  @test boundingbox(m) == Box(point(-1, -1), point(1, 1))
  @test boundingbox(d) == Box(point(-1, -1), point(1, 1))
  @test @allocated(boundingbox(m)) < 50
  @test @allocated(boundingbox(d)) < 50

  d = PointSet(T[0 1 2; 0 2 1])
  @test boundingbox(d) == Box(point(0, 0), point(2, 2))
  @test @allocated(boundingbox(d)) < 50
  d = PointSet(T[1 2; 2 1])
  @test boundingbox(d) == Box(point(1, 1), point(2, 2))
  @test @allocated(boundingbox(d)) < 50

  d = cartgrid(10, 10)
  @test boundingbox(d) == Box(point(0, 0), point(10, 10))
  @test @allocated(boundingbox(d)) < 50
  d = cartgrid(100, 200)
  @test boundingbox(d) == Box(point(0, 0), point(100, 200))
  @test @allocated(boundingbox(d)) < 50
  d = CartesianGrid((10, 10), T.((1, 1)), T.((1, 1)))
  @test boundingbox(d) == Box(point(1, 1), point(11, 11))
  @test @allocated(boundingbox(d)) < 50

  d = PointSet(T[0 1 2; 0 2 1])
  v = view(d, 1:2)
  @test boundingbox(v) == Box(point(0, 0), point(1, 2))
  @test @allocated(boundingbox(v)) < 50

  d = cartgrid(10, 10)
  v = view(d, 1:2)
  @test boundingbox(v) == Box(point(0, 0), point(2, 1))
  @test @allocated(boundingbox(v)) < 9000

  g = cartgrid(10, 10)
  d = convert(RectilinearGrid, g)
  @test boundingbox(d) == Box(point(0, 0), point(10, 10))
  @test @allocated(boundingbox(d)) < 50

  g = cartgrid(10, 10)
  d = TransformedGrid(g, Rotate(T(π / 2)))
  @test boundingbox(d) ≈ Box(point(-10, 0), point(0, 10))
  @test @allocated(boundingbox(d)) < 3000

  g = cartgrid(10, 10)
  rg = convert(RectilinearGrid, g)
  d = TransformedGrid(rg, Rotate(T(π / 2)))
  @test boundingbox(d) ≈ Box(point(-10, 0), point(0, 10))
  @test @allocated(boundingbox(d)) < 3000

  g = cartgrid(10, 10)
  m = convert(SimpleMesh, g)
  @test boundingbox(m) == Box(point(0, 0), point(10, 10))
  @test @allocated(boundingbox(m)) < 50

  p = ParaboloidSurface(point(1, 2, 3), T(5), T(4))
  @test boundingbox(p) ≈ Box(point(-4, -3, 3), point(6, 7, 73 / 16))
end
