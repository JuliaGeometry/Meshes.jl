@testset "orientation" begin
  # test orientation
  t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
  @test orientation(t) == CCW
  t = Triangle(P2(0, 0), P2(0, 1), P2(1, 0))
  @test orientation(t) == CW

  # orientation of 3D rings in X-Y plane
  r1 = Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  r2 = Ring(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
  @test orientation(r1) == orientation(r2)
end
