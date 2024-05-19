@testset "sideof" begin
  p1, p2, p3 = point(0, 0), point(1, 1), point(0.25, 0.5)
  l = Line(point(0.5, 0.0), point(0.0, 1.0))
  @test sideof(p1, l) == LEFT
  @test sideof(p2, l) == RIGHT
  @test sideof(p3, l) == ON
  pts = [p1, p2, p3]
  @test sideof(pts, l) == [LEFT, RIGHT, ON]

  p1, p2, p3 = point(0.5, 0.5), point(1.5, 0.5), point(1, 1)
  c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test sideof(p1, c) == IN
  @test sideof(p2, c) == OUT
  @test sideof(p3, c) == IN
  pts = [p1, p2, p3]
  @test sideof(pts, c) == [IN, OUT, IN]

  points = point.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1)])
  connec = connect.([(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)], Triangle)
  mesh = SimpleMesh(points, connec)
  @test sideof(point(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(point(0.25, 0.25, -0.1), mesh) == OUT
  pts = point.([(0.25, 0.25, 0.1), (0.25, 0.25, -0.1)])
  @test sideof(pts, mesh) == [IN, OUT]

  # ray goes through vertex
  @test sideof(point(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(point(0.25, 0.25, -0.1), mesh) == OUT

  # ray goes through edge of triangle
  @test sideof(point(0.1, 0.1, 0.1), mesh) == IN
  @test sideof(point(0.1, 0.1, -0.1), mesh) == OUT

  # point coincides with edge of triangle
  @test sideof(point(0.5, 0.0, 0.0), mesh) == IN

  # point coincides with corner of triangle
  @test sideof(point(0.0, 0.0, 0.0), mesh) == IN

  points = point.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])
  mesh = SimpleMesh(points, connec)
  # ray collinear with edge
  @test sideof(point(0.0, 0.0, 0.1), mesh) == IN
  @test sideof(point(0.0, 0.0, -0.1), mesh) == OUT

  # sideof for meshes that have elements > 3-gons.
  points = point.([(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1), (1, 1, 0)])
  connec = connect.([(1, 2, 4), (1, 4, 3), (2, 3, 4), (1, 2, 5, 3)])
  mesh = SimpleMesh(points, connec)
  @test sideof(point(0.25, 0.25, 0.1), mesh) == IN

  # sideof only defined for surface meshes
  points = point.([(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)])
  connec = connect.([(1, 2, 3, 4)], [Tetrahedron])
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError("winding number only defined for surface meshes") sideof(point(0, 0, 0), mesh)
end
