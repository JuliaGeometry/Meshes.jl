@testset "Distances" begin
  p = point(0, 1)
  l = Line(point(0, 0), point(1, 0))
  @test evaluate(Euclidean(), p, l) == T(1) * u"m"
  @test evaluate(Euclidean(), l, p) == T(1) * u"m"

  p = point(68, 259)
  l = Line(point(68, 260), point(69, 261))
  @test evaluate(Euclidean(), p, l) ≤ T(0.8) * u"m"

  line1 = Line(point(-1, 0, 0), point(1, 0, 0))
  line2 = Line(point(0, -1, 1), point(0, 1, 1))  # line2 ⟂ line1, z++
  line3 = Line(point(-1, 1, 0), point(1, 1, 0))  # line3 ∥ line1
  line4 = Line(point(-2, 0, 0), point(2, 0, 0))  # line4 colinear with line1
  line5 = Line(point(0, -1, 0), point(0, 1, 0))  # line5 intersects line1
  @test evaluate(Euclidean(), line1, line2) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line3) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line4) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), line1, line5) ≈ T(0) * u"m"

  p1, p2 = point(1, 0), point(0, 1)
  @test evaluate(Chebyshev(), p1, p2) == T(1) * u"m"
  @test evaluate(Euclidean(), p1, p2) == T(√2) * u"m"

  latlon1 = LatLon(T(0), T(0))
  latlon2 = LatLon(T(1), T(0))
  cart1 = convert(Cartesian, latlon1)
  cart2 = convert(Cartesian, latlon2)
  p1 = Point(latlon1)
  p2 = Point(latlon2)
  p3 = Point(cart1)
  p4 = Point(cart2)
  @test evaluate(Haversine(), p1, p2) ≈ T(111194.92664455874) * u"m"
  @test evaluate(Haversine(), p3, p4) ≈ T(111194.92664455874) * u"m"
  @test evaluate(SphericalAngle(), p1, p2) ≈ deg2rad(T(1) * u"°")
  @test evaluate(SphericalAngle(), p3, p4) ≈ deg2rad(T(1) * u"°")
end
