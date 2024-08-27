@testitem "orientation" setup = [Setup] begin
  # test orientation
  t = Triangle(cart(0, 0), cart(1, 0), cart(0, 1))
  @test orientation(t) == CCW
  t = Triangle(cart(0, 0), cart(0, 1), cart(1, 0))
  @test orientation(t) == CW

  # orientation of 3D rings in X-Y plane
  r1 = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  r2 = Ring(cart.([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]))
  @test orientation(r1) == orientation(r2)

  # orientation of rings with LatLon coordinates
  r = Ring(latlon.([(0, 0), (0, 90), (90, 0)]))
  @test orientation(r) == CCW
  r = Ring(latlon.([(0, 0), (90, 0), (0, 90)]))
  @test orientation(r) == CW
end
