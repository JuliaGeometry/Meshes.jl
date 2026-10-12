@testitem "Booleans" setup = [Setup] begin
  # union of overlapping polygons
  q1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  q2 = Quadrangle(cart(2, 2), cart(6, 2), cart(6, 6), cart(2, 6))
  u = q1 ∪ q2
  @test measure(u) ≈ T(28) * u"m^2"
  @test all(
    vertices(u) .≈ [cart(2, 4), cart(0, 4), cart(0, 0), cart(4, 0), cart(4, 2), cart(6, 2), cart(6, 6), cart(2, 6)]
  )

  # difference of overlapping polygons
  q1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  q2 = Quadrangle(cart(2, 2), cart(6, 2), cart(6, 6), cart(2, 6))
  d1 = setdiff(q1, q2)
  d2 = setdiff(q2, q1)
  @test measure(d1) ≈ T(12) * u"m^2"
  @test measure(d2) ≈ T(12) * u"m^2"
  @test all(vertices(d1) .≈ [cart(2, 4), cart(0, 4), cart(0, 0), cart(4, 0), cart(4, 2), cart(2, 2)])

  # non-convex polygons
  p1 = Pgon(cart.([(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)]))
  p2 = Pgon(cart.([(0, 0), (4, 0), (4, 4), (3, 4), (3, 1), (0, 1)]))
  @test measure(p1 ∪ p2) ≈ T(10) * u"m^2"
  @test measure(setdiff(p1, p2)) ≈ T(3) * u"m^2"

  # polygon inside the other one
  q1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  q2 = Quadrangle(cart(1, 1), cart(3, 1), cart(3, 3), cart(1, 3))
  u = q1 ∪ q2
  d = setdiff(q1, q2)
  @test measure(u) ≈ T(16) * u"m^2"
  @test length(rings(u)) == 1
  @test measure(d) ≈ T(12) * u"m^2"
  @test length(rings(d)) == 2

  # disjoint polygons
  q1 = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  q2 = Quadrangle(cart(10, 10), cart(11, 10), cart(11, 11), cart(10, 11))
  u = q1 ∪ q2
  d = setdiff(q1, q2)
  @test u isa Multi
  @test measure(u) ≈ T(17) * u"m^2"
  @test measure(d) ≈ T(16) * u"m^2"

  # difference with itself
  q = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  @test isnothing(setdiff(q, q))

  # CRS propagation
  q1 = Quadrangle(merc(0, 0), merc(4, 0), merc(4, 4), merc(0, 4))
  q2 = Quadrangle(merc(2, 2), merc(6, 2), merc(6, 6), merc(2, 6))
  @test crs(q1 ∪ q2) === crs(q1)
  @test crs(setdiff(q1, q2)) === crs(q1)
end
