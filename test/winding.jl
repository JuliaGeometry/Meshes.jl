@testset "winding" begin
  p = point(0.5, 0.5)
  c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]

  p = point(0.5, 0.5)
  c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)]))
  @test winding(p, c) ≈ T(2)
  @test winding(p, reverse(c)) ≈ T(-2)
  @test winding([p, p], c) ≈ T[2, 2]

  p = Point(LatLon(T(0.5), T(0.5)))
  c = Ring(Point.([LatLon(T(0), T(0)), LatLon(T(1), T(0)), LatLon(T(1), T(1)), LatLon(T(0), T(1))]))
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]

  # error: both arguments must have the same CRS
  p = Point(LatLon(T(0.5), T(0.5)))
  c = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)]))
  @test_throws ArgumentError winding(p, c)

  m = boundary(Box(point(0, 0, 0), point(2, 2, 2)))
  @test all(>(0), winding(vertices(m), m))
  @test isapprox(winding(point(1, 1, 1), m), T(1), atol=atol(T))
  @test isapprox(winding(point(3, 3, 3), m), T(0), atol=atol(T))
end
