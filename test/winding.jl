@testset "winding" begin
  p = P2(0.5, 0.5)
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  @test winding(p, c) ≈ T(1)
  @test winding(p, reverse(c)) ≈ T(-1)
  @test winding([p, p], c) ≈ T[1, 1]

  p = P2(0.5, 0.5)
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)])
  @test winding(p, c) ≈ T(2)
  @test winding(p, reverse(c)) ≈ T(-2)
  @test winding([p, p], c) ≈ T[2, 2]

  m = boundary(Box(P3(0, 0, 0), P3(2, 2, 2)))
  @test all(>(0), winding(vertices(m), m))
  @test isapprox(winding(P3(1, 1, 1), m), T(1), atol = atol(T))
  @test isapprox(winding(P3(3, 3, 3), m), T(0), atol = atol(T))
end
