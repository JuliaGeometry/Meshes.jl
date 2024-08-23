@testitem "Bounding boxes" begin
  p = cart(0, 0)
  @test boundingbox(p) == Box(p, p)
  @test @allocated(boundingbox(p)) < 50

  b = Box(cart(0, 0), cart(1, 1))
  @test boundingbox(b) == b
  @test @allocated(boundingbox(b)) < 50

  r = Ray(cart(0, 0), vector(1, 0))
  @test boundingbox(r) == Box(cart(0, 0), cart(T(Inf), 0))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(cart(1, 1), vector(0, 1))
  @test boundingbox(r) == Box(cart(1, 1), cart(1, T(Inf)))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(cart(1, 1), vector(-1, -1))
  @test boundingbox(r) == Box(cart(T(-Inf), T(-Inf)), cart(1, 1))
  @test @allocated(boundingbox(r)) < 50
  r = Ray(cart(-1, 1), vector(1, -1))
  @test boundingbox(r) == Box(cart(-1, T(-Inf)), cart(T(Inf), 1))
  @test @allocated(boundingbox(r)) < 50

  b = Ball(cart(0, 0), T(1))
  @test boundingbox(b) == Box(cart(-1, -1), cart(1, 1))
  @test @allocated(boundingbox(b)) < 50
  b = Ball(cart(1, 1), T(1))
  @test boundingbox(b) == Box(cart(0, 0), cart(2, 2))
  @test @allocated(boundingbox(b)) < 50

  s = Sphere(cart(0, 0), T(1))
  @test boundingbox(s) == Box(cart(-1, -1), cart(1, 1))
  @test @allocated(boundingbox(s)) < 50
  s = Sphere(cart(1, 1), T(1))
  @test boundingbox(s) == Box(cart(0, 0), cart(2, 2))
  @test @allocated(boundingbox(s)) < 50

  c = Cylinder(T(1))
  b = boundingbox(c)
  @test b == Box(cart(-1, -1, 0), cart(1, 1, 1))

  c = CylinderSurface(T(1))
  b = boundingbox(c)
  @test b == Box(cart(-1, -1, 0), cart(1, 1, 1))

  c = Cone(Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(1)), cart(0, 0, 1))
  b = boundingbox(c)
  @test b == Box(cart(-1, -1, 0), cart(1, 1, 1))

  c = ConeSurface(Disk(Plane(cart(0, 0, 0), vector(0, 0, 1)), T(1)), cart(0, 0, 1))
  b = boundingbox(c)
  @test b == Box(cart(-1, -1, 0), cart(1, 1, 1))

  b = Box(cart(-3, -1), cart(0.5, 0.5))
  s = Sphere(cart(0, 0), T(2))
  m = Multi([b, s])
  d = GeometrySet([b, s])
  @test boundingbox(m) == Box(cart(-3, -2), cart(2, 2))
  @test boundingbox(d) == Box(cart(-3, -2), cart(2, 2))
  @test @allocated(boundingbox(m)) < 2500
  @test @allocated(boundingbox(d)) < 2500

  b1 = Box(cart(0, 0), cart(1, 1))
  b2 = Box(cart(-1, -1), cart(0.5, 0.5))
  m = Multi([b1, b2])
  d = GeometrySet([b1, b2])
  @test boundingbox(m) == Box(cart(-1, -1), cart(1, 1))
  @test boundingbox(d) == Box(cart(-1, -1), cart(1, 1))
  @test @allocated(boundingbox(m)) < 50
  @test @allocated(boundingbox(d)) < 50

  d = PointSet(cart(0, 0), cart(1, 2), cart(2, 1))
  @test boundingbox(d) == Box(cart(0, 0), cart(2, 2))
  @test @allocated(boundingbox(d)) < 50
  d = PointSet(cart(1, 2), cart(2, 1))
  @test boundingbox(d) == Box(cart(1, 1), cart(2, 2))
  @test @allocated(boundingbox(d)) < 50

  d = cartgrid(10, 10)
  @test boundingbox(d) == Box(cart(0, 0), cart(10, 10))
  @test @allocated(boundingbox(d)) < 50
  d = cartgrid(100, 200)
  @test boundingbox(d) == Box(cart(0, 0), cart(100, 200))
  @test @allocated(boundingbox(d)) < 50
  d = CartesianGrid((10, 10), T.((1, 1)), T.((1, 1)))
  @test boundingbox(d) == Box(cart(1, 1), cart(11, 11))
  @test @allocated(boundingbox(d)) < 50

  d = PointSet(cart(0, 0), cart(1, 2), cart(2, 1))
  v = view(d, 1:2)
  @test boundingbox(v) == Box(cart(0, 0), cart(1, 2))
  @test @allocated(boundingbox(v)) < 50

  d = cartgrid(10, 10)
  v = view(d, 1:2)
  @test boundingbox(v) == Box(cart(0, 0), cart(2, 1))
  @test @allocated(boundingbox(v)) < 10000

  g = cartgrid(10, 10)
  d = convert(RectilinearGrid, g)
  @test boundingbox(d) == Box(cart(0, 0), cart(10, 10))
  @test @allocated(boundingbox(d)) < 50

  g = cartgrid(10, 10)
  d = TransformedGrid(g, Rotate(T(π / 2)))
  @test boundingbox(d) ≈ Box(cart(-10, 0), cart(0, 10))
  @test @allocated(boundingbox(d)) < 3000

  g = cartgrid(10, 10)
  rg = convert(RectilinearGrid, g)
  d = TransformedGrid(rg, Rotate(T(π / 2)))
  @test boundingbox(d) ≈ Box(cart(-10, 0), cart(0, 10))
  @test @allocated(boundingbox(d)) < 3000

  g = cartgrid(10, 10)
  m = convert(SimpleMesh, g)
  @test boundingbox(m) == Box(cart(0, 0), cart(10, 10))
  @test @allocated(boundingbox(m)) < 50

  p = ParaboloidSurface(cart(1, 2, 3), T(5), T(4))
  @test boundingbox(p) ≈ Box(cart(-4, -3, 3), cart(6, 7, 73 / 16))

  # CRS propagation
  r = Ray(merc(-1, 1), vector(1, -1))
  @test crs(boundingbox(r)) === crs(r)
  g = CartesianGrid((10, 10), merc(0, 0), (T(1), T(1)))
  m = convert(SimpleMesh, g)
  @test crs(boundingbox(m)) === crs(m)
end
