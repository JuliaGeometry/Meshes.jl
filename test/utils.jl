@testset "Utilities" begin
  a, b, c = point(0, 0), point(1, 0), point(0, 1)
  @test signarea(a, b, c) == T(0.5) * u"m^2"
  a, b, c = point(0, 0), point(0, 1), point(1, 0)
  @test signarea(a, b, c) == T(-0.5) * u"m^2"

  normals = [
    vector(1, 0, 0),
    vector(0, 1, 0),
    vector(0, 0, 1),
    vector(-1, 0, 0),
    vector(0, -1, 0),
    vector(0, 0, -1),
    vector(ntuple(i -> rand() - 0.5, 3))
  ]
  for n in normals
    u, v = householderbasis(n)
    @test u isa Vec{3}
    @test v isa Vec{3}
    @test ustrip.(u × v) ≈ n ./ norm(n)
  end
  n = Vec(T(1) * u"cm", T(1) * u"cm", T(1) * u"cm")
  u, v = householderbasis(n)
  @test unit(eltype(u)) == u"cm"
  @test unit(eltype(v)) == u"cm"
  n = Vec(T(1) * u"km", T(1) * u"km", T(1) * u"km")
  u, v = householderbasis(n)
  @test unit(eltype(u)) == u"km"
  @test unit(eltype(v)) == u"km"

  @test Meshes.mayberound(1.1, 1.0, 0.2) ≈ 1.0
  @test Meshes.mayberound(1.1, 1.0, 0.10000000000000001) ≈ 1.1
  @test Meshes.mayberound(1.1, 1.0, 0.05) ≈ 1.1

  # intersect parameters
  p1, p2 = point(0, 0), point(1, 1)
  p3, p4 = point(1, 0), point(0, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1, p2 = point(0, 0, 0), point(1, 1, 1)
  p3, p4 = point(1, 0, 0), point(0, 1, 1)
  @inferred Meshes.intersectparameters(p1, p2, p3, p4)
  @inferred Meshes.intersectparameters(p1, p3, p2, p4)
  @inferred Meshes.intersectparameters(p1, p2, p1, p2)

  p1 = point(1, 1)
  p2 = Point(Polar(T(√2), T(π / 4)))
  p3 = Point(LatLon(T(30), T(60)))
  p4 = Point(GeocentricLatLon(T(30), T(60)))
  p5 = Point(AuthalicLatLon(T(30), T(60)))
  p6 = Point(Mercator(T(1), T(1)))
  c1 = Rope(Point(Polar(T(0), T(0))), Point(Polar(T(1), T(0))), Point(Polar(T(1), T(π / 2))))
  c2 = Ring(Point(Polar(T(0), T(0))), Point(Polar(T(1), T(0))), Point(Polar(T(1), T(π / 2))))
  t = Triangle(Point(Polar(T(0), T(0))), Point(Polar(T(1), T(0))), Point(Polar(T(1), T(π / 2))))
  @test Meshes.ascart2(p1) == point(1, 1)
  @test Meshes.ascart2(p2) ≈ point(1, 1)
  @test Meshes.ascart2(p3) == Point(Cartesian{WGS84Latest}(T(30), T(60)))
  @test Meshes.ascart2(p5) == Point(Cartesian{WGS84Latest}(T(30), T(60)))
  @test Meshes.ascart2(p5) == Point(Cartesian{WGS84Latest}(T(30), T(60)))
  @test Meshes.ascart2(p6) == Point(Cartesian{WGS84Latest}(T(1), T(1)))
  @test Meshes.ascart2(c1) ≈ Rope(point(0, 0), point(1, 0), point(0, 1))
  @test Meshes.ascart2(c2) ≈ Ring(point(0, 0), point(1, 0), point(0, 1))
  @test Meshes.ascart2(t) ≈ Triangle(point(0, 0), point(1, 0), point(0, 1))
  # error: points must have 2 coordinates
  @test_throws ArgumentError Meshes.ascart2(point(1, 1, 1))
end
