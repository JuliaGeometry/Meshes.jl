@testset "orientation" begin
  # test orientation
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  @test orientation(t) == CCW
  t = Triangle(cart(0, 0), cart(0, 1), cart(1, 0))
  @test orientation(t) == CW

  # orientation of 3D rings in X-Y plane
  r1 = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  r2 = Ring(cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
  @test orientation(r1) == orientation(r2)
end
