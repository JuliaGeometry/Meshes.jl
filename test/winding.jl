@testset "winding" begin
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  @test winding(P2(0.5, 0.5), c) ≈ 1
  @test winding(P2(0.5, 0.5), reverse(c)) ≈ -1
  c = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0), (1, 0), (1, 1), (0, 1)])
  @test winding(P2(0.5, 0.5), c) ≈ 2
  @test winding(P2(0.5, 0.5), reverse(c)) ≈ -2
end
