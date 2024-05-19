@testset "orientation" begin
  # test orientation
  t = Triangle(point(0, 0), point(1, 0), point(0, 1))
  @test orientation(t) == CCW
  t = Triangle(point(0, 0), point(0, 1), point(1, 0))
  @test orientation(t) == CW

  # orientation of 3D rings in X-Y plane
  r1 = Ring(point.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  r2 = Ring(point.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
  @test orientation(r1) == orientation(r2)
end
