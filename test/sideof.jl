@testset "sideof" begin
  p1, p2, p3 = cart(0, 0), cart(1, 1), cart(0.25, 0.5)
  l = Line(cart(0.5, 0.0), cart(0.0, 1.0))
  @test sideof(p1, l) == LEFT
  @test sideof(p2, l) == RIGHT
  @test sideof(p3, l) == ON
  pts = [p1, p2, p3]
  @test sideof(pts, l) == [LEFT, RIGHT, ON]

  p1, p2, p3 = cart(0.5, 0.5), cart(1.5, 0.5), cart(1, 1)
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test sideof(p1, c) == IN
  @test sideof(p2, c) == OUT
  @test sideof(p3, c) == ON
  pts = [p1, p2, p3]
  @test sideof(pts, c) == [IN, OUT, ON]

  c_mm = LengthUnit(Unitful.mm)(c)
  @test sideof(p1, c_mm) == IN
  @test sideof(p2, c_mm) == OUT
  @test sideof(p3, c_mm) == ON
  pts = [p1, p2, p3]
  @test sideof(pts, c_mm) == [IN, OUT, ON]

  p1, p2, p3 = merc(0.5, 0.5), merc(1.5, 0.5), merc(1, 1)
  c = Ring([merc(0, 0), merc(1, 0), merc(1, 1), merc(0, 1)])
  @test sideof(p1, c) == IN
  @test sideof(p2, c) == OUT
  @test sideof(p3, c) == ON
  pts = [p1, p2, p3]
  @test sideof(pts, c) == [IN, OUT, ON]

  points = cart.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1)])
  connec = connect.([(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)], Triangle)
  mesh = SimpleMesh(points, connec)
  @test sideof(cart(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(cart(0.25, 0.25, -0.1), mesh) == OUT
  pts = cart.([(0.25, 0.25, 0.1), (0.25, 0.25, -0.1)])
  @test sideof(pts, mesh) == [IN, OUT]

  # ray goes through vertex
  @test sideof(cart(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(cart(0.25, 0.25, -0.1), mesh) == OUT

  # ray goes through edge of triangle
  @test sideof(cart(0.1, 0.1, 0.1), mesh) == IN
  @test sideof(cart(0.1, 0.1, -0.1), mesh) == OUT

  # point coincides with edge of triangle
  @test sideof(cart(0.5, 0.0, 0.0), mesh) == IN

  # point coincides with corner of triangle
  @test sideof(cart(0.0, 0.0, 0.0), mesh) == IN

  points = cart.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
  mesh = SimpleMesh(points, connec)
  # ray collinear with edge
  @test sideof(cart(0.0, 0.0, 0.1), mesh) == IN
  @test sideof(cart(0.0, 0.0, -0.1), mesh) == OUT

  # sideof for meshes that have elements > 3-gons.
  points = cart.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1), (1, 1, 0)])
  connec = connect.([(1, 2, 4), (1, 4, 3), (2, 3, 4), (1, 2, 5, 3)])
  mesh = SimpleMesh(points, connec)
  @test sideof(cart(0.25, 0.25, 0.1), mesh) == IN

  # sideof only defined for surface meshes
  points = cart.([(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)])
  connec = connect.([(1, 2, 3, 4)], [Tetrahedron])
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError("winding number only defined for surface meshes") sideof(cart(0, 0, 0), mesh)

  # sideof serial vs threads
  p = first(randpoint2(1))
  r = Ring(randpoint2(500))
  serial = Meshes._sideofserial(p, r)
  threads = Meshes._sideofthreads(p, r)
  @test serial == threads
end
