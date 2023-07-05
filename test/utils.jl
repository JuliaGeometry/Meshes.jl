@testset "Utilities" begin
  a, b, c = P2(0, 0), P2(1, 0), P2(0, 1)
  @test signarea(a, b, c) == T(0.5)
  a, b, c = P2(0, 0), P2(0, 1), P2(1, 0)
  @test signarea(a, b, c) == T(-0.5)

  p1, p2, p3 = P2(0, 0), P2(1, 1), P2(0.25, 0.5)
  s = Segment(P2(0.5, 0.0), P2(0.0, 1.0))
  @test sideof(p1, s) == :LEFT
  @test sideof(p2, s) == :RIGHT
  @test sideof(p3, s) == :ON

  p1, p2, p3 = P2(0.5, 0.5), P2(1.5, 0.5), P2(1, 1)
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  @test sideof(p1, c) == :INSIDE
  @test sideof(p2, c) == :OUTSIDE
  @test sideof(p3, c) == :INSIDE

  @test iscollinear(P2(0, 0), P2(1, 1), P2(2, 2))
  @test iscoplanar(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0))

  # drop units from unitful value and type
  @test Meshes.dropunits(1.0u"mm") == Float64
  @test Meshes.dropunits(typeof(1.0u"mm")) == Float64

  # return the same type in case of no units
  @test Meshes.dropunits(1.0) == Float64
  @test Meshes.dropunits(Float64) == Float64

  normals = [V3(1, 0, 0), V3(0, 1, 0), V3(0, 0, 1), V3(-1, 0, 0), V3(0, -1, 0), V3(0, 0, -1), V3(rand(3) .- 0.5)]
  for n in normals
    u, v = householderbasis(n)
    @test u × v ≈ n ./ norm(n)
  end

  @test mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test mayberound(1.1, 1.0, 0.05) ≈ 1.1

  # point in mesh
  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1)]
  connec = connect.([(1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4)], Triangle)
  mesh = SimpleMesh(points, connec)
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == :INSIDE
  @test sideof(P3(0.25, 0.25, -0.1), mesh) == :OUTSIDE

  # ray goes through vertex
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == :INSIDE
  @test sideof(P3(0.25, 0.25, -0.1), mesh) == :OUTSIDE

  # ray goes through edge of triangle
  @test sideof(P3(0.1, 0.1, 0.1), mesh) == :INSIDE
  @test sideof(P3(0.1, 0.1, -0.1), mesh) == :OUTSIDE

  # point coincides with edge of triangle
  @test sideof(P3(0.5, 0.0, 0.0), mesh) == :ON

  # point coincides with corner of triangle
  @test sideof(P3(0.0, 0.0, 0.0), mesh) == :ON

  # point on face of triangle
  @test sideof(P3(0.1, 0.1, 0.0), mesh) == :ON

  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
  mesh = SimpleMesh(points, connec)
  # ray collinear with edge
  @test sideof(P3(0.0, 0.0, 0.1), mesh) == :INSIDE
  @test sideof(P3(0.0, 0.0, -0.1), mesh) == :OUTSIDE

  # sideof for meshes that have elements > 3-gons.
  points = P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0.25, 0.25, 1), (1, 1, 0)]
  connec = connect.([(1, 2, 4), (1, 4, 3), (2, 3, 4), (1, 2, 5, 3)])
  mesh = SimpleMesh(points, connec)
  @test sideof(P3(0.25, 0.25, 0.1), mesh) == :INSIDE

  # sideof only defined for surface meshes
  points = P3[(0, 0, 0), (1, 0, 0), (1, 1, 1), (0, 1, 0)]
  connec = connect.([(1, 2, 3, 4)], [Tetrahedron])
  mesh = SimpleMesh(points, connec)
  @test_throws AssertionError("sideof only defined for surface meshes") sideof(P3(0, 0, 0), mesh)

  # intersect parameters
  p1, p2 = P2(0, 0), P2(1, 1)
  p3, p4 = P2(1, 0), P2(0, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1, p2 = P3(0, 0, 0), P3(1, 1, 1)
  p3, p4 = P3(1, 0, 0), P3(0, 1, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  # overlap parameters
  @test Meshes.overlapparameters(2.5, 1.4, 1.1, 2.0) == (1.4, 2.0)
  @test Meshes.overlapparameters(2.0, 1.1, 1.4, 2.5) == (1.4, 2.0)
  @test Meshes.overlapparameters(2.0, 2.5, 1.1, 1.4) == (1.4, 2.0)
end
