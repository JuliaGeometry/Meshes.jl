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
end
