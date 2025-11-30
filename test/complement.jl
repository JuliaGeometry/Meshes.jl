@testitem "Complement of geometries" setup = [Setup] begin
  τ = Meshes.atol(T)

  t = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  p = !t
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(cart.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(cart.([(0, 0), (1, 1), (1, 0)]))

  q = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  p = !q
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(cart.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  p = PolyArea(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  n = !p
  r = rings(n)
  @test n isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(cart.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  o = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  i1 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
  i2 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
  p = PolyArea([o, reverse(i1), reverse(i2)])
  m = !p
  r = rings(m)
  @test m isa MultiPolygon
  @test length(r) == 4
  g = parent(m)
  @test length(g) == 3
  @test rings(g[1])[1] ≈ Ring(cart.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test rings(g[1])[2] == Ring(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))
  @test rings(g[2]) == [Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))]
  @test rings(g[3]) == [Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))]

  b = Box(cart(0, 0), cart(1, 1))
  p = !b
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(cart.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(cart.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  b = Ball(cart(0, 0), T(1))
  p = !b
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(cart.([(-1 - τ, -1 - τ), (1 + τ, -1 - τ), (1 + τ, 1 + τ), (-1 - τ, 1 + τ)]))
end
