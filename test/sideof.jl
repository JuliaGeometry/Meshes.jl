@testset "sideof" begin
  # -----
  # LINE
  # -----

  p1, p2, p3 = P2(0, 0), P2(1, 1), P2(0.25, 0.5)
  l = Line(P2(0.5, 0.0), P2(0.0, 1.0))
  @test sideof(p1, l) == LEFT
  @test sideof(p2, l) == RIGHT
  @test sideof(p3, l) == ON

  # -----
  # RING
  # -----

  p1, p2, p3 = P2(0.5, 0.5), P2(1.5, 0.5), P2(1, 1)
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  @test sideof(p1, c) == IN
  @test sideof(p2, c) == OUT
  @test sideof(p3, c) == IN

  # -----
  # MESH
  # -----

  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1)]
  connec = connect.([(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)], Triangle)
  mesh = SimpleMesh(points, connec)
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(P3(0.25, 0.25, -0.1), mesh) == OUT

  # ray goes through vertex
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == IN
  @test sideof(P3(0.25, 0.25, -0.1), mesh) == OUT

  # ray goes through edge of triangle
  @test sideof(P3(0.1, 0.1, 0.1), mesh) == IN
  @test sideof(P3(0.1, 0.1, -0.1), mesh) == OUT

  # point coincides with edge of triangle
  @test sideof(P3(0.5, 0.0, 0.0), mesh) == ON

  # point coincides with corner of triangle
  @test sideof(P3(0.0, 0.0, 0.0), mesh) == ON

  # point on face of triangle
  @test sideof(P3(0.1, 0.1, 0.0), mesh) == ON

  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
  mesh = SimpleMesh(points, connec)
  # ray collinear with edge
  @test sideof(P3(0.0, 0.0, 0.1), mesh) == IN
  @test sideof(P3(0.0, 0.0, -0.1), mesh) == OUT

  # sideof for meshes that have elements > 3-gons.
  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1), (1, 1, 0)]
  connec = connect.([(1, 2, 4), (1, 4, 3), (2, 3, 4), (1, 2, 5, 3)])
  mesh = SimpleMesh(points, connec)
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == IN

  # sideof only defined for surface meshes
  points = P3[(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)]
  connec = connect.([(1, 2, 3, 4)], Tetrahedron)
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError("sideof only defined for surface meshes") sideof(P3(0, 0, 0), mesh)
end
