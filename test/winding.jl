@testitem "Winding numbers" setup = [Setup] begin
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

  p = merc(0.5, 0.5)
  c = Ring([merc(0, 0), merc(1, 0), merc(1, 1), merc(0, 1)])
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]
  # record allocations for merc
  allocmerc = @allocated winding(p, c)

  # exact same memory allocations
  @test alloccart == allocmerc

  m = boundary(Box(cart(0, 0, 0), cart(2, 2, 2)))
  @test all(>(0), winding(vertices(m), m))
  @test isapprox(winding(cart(1, 1, 1), m), T(1), atol=Meshes.atol(T))
  @test isapprox(winding(cart(3, 3, 3), m), T(0), atol=Meshes.atol(T))
end
