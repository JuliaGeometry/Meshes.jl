@testset "winding" begin
  p = cart(0.5, 0.5)
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]

  p = cart(0.5, 0.5)
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)]))
  @test winding(p, c) ≈ T(2)
  @test winding(p, reverse(c)) ≈ T(-2)
  @test winding([p, p], c) ≈ T[2, 2]
  # record allocations for cartesian
  alloccart = @allocated winding(p, c)

  p = latlon(0.5, 0.5)
  c = Ring([latlon(0, 0), latlon(1, 0), latlon(1, 1), latlon(0, 1)])
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]
  # record allocations for latlon
  alloclatlon = @allocated winding(p, c)

  # exact same memory allocations
  @test alloccart == alloclatlon

  m = boundary(Box(cart(0, 0, 0), cart(2, 2, 2)))
  @test all(>(0), winding(vertices(m), m))
  @test isapprox(winding(cart(1, 1, 1), m), T(1), atol=atol(T))
  @test isapprox(winding(cart(3, 3, 3), m), T(0), atol=atol(T))
end
