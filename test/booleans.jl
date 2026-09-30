@testitem "Booleans" setup = [Setup] begin
  sq₁ = Quadrangle(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  sq₂ = Quadrangle(cart(2, 2), cart(6, 2), cart(6, 6), cart(2, 6))

  # union of overlapping polygons
  u = union(sq₁, sq₂)
  @test measure(u) ≈ T(28) * u"m^2"
  @test all(
    vertices(u) .≈ [cart(2, 4), cart(0, 4), cart(0, 0), cart(4, 0), cart(4, 2), cart(6, 2), cart(6, 6), cart(2, 6)]
  )

  # difference of overlapping polygons
  d = setdiff(sq₁, sq₂)
  @test measure(d) ≈ T(12) * u"m^2"
  @test all(vertices(d) .≈ [cart(2, 4), cart(0, 4), cart(0, 0), cart(4, 0), cart(4, 2), cart(2, 2)])
  @test measure(setdiff(sq₂, sq₁)) ≈ T(12) * u"m^2"

  # non-convex polygons
  poly₁ = PolyArea(cart.([(0, 0), (4, 0), (4, 1), (1, 1), (1, 4), (0, 4)]))
  poly₂ = PolyArea(cart.([(0, 0), (4, 0), (4, 4), (3, 4), (3, 1), (0, 1)]))
  @test measure(union(poly₁, poly₂)) ≈ T(10) * u"m^2"
  @test measure(setdiff(poly₁, poly₂)) ≈ T(3) * u"m^2"

  # polygon inside the other one
  inner = Quadrangle(cart(1, 1), cart(3, 1), cart(3, 3), cart(1, 3))
  @test measure(union(sq₁, inner)) ≈ T(16) * u"m^2"
  @test length(rings(union(sq₁, inner))) == 1
  d = setdiff(sq₁, inner)
  @test measure(d) ≈ T(12) * u"m^2"
  @test length(rings(d)) == 2

  # disjoint polygons
  far = Quadrangle(cart(10, 10), cart(11, 10), cart(11, 11), cart(10, 11))
  u = union(sq₁, far)
  @test u isa Multi
  @test measure(u) ≈ T(17) * u"m^2"
  @test measure(setdiff(sq₁, far)) ≈ T(16) * u"m^2"

  # difference with itself
  @test isnothing(setdiff(sq₁, sq₁))

  # CRS propagation
  poly₁ = Quadrangle(merc(0, 0), merc(4, 0), merc(4, 4), merc(0, 4))
  poly₂ = Quadrangle(merc(2, 2), merc(6, 2), merc(6, 6), merc(2, 6))
  @test crs(union(poly₁, poly₂)) === crs(poly₁)
  @test crs(setdiff(poly₁, poly₂)) === crs(poly₁)
end
