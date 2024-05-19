@testset "complement" begin
  τ = atol(T)

  t = Triangle(point(0, 0), point(1, 0), point(1, 1))
  p = !t
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(point.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(point.([(0, 0), (1, 1), (1, 0)]))

  q = Quadrangle(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
  p = !q
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(point.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  p = PolyArea(point(0, 0), point(1, 0), point(1, 1), point(0, 1))
  n = !p
  r = rings(n)
  @test n isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(point.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  o = point.([(0, 0), (1, 0), (1, 1), (0, 1)])
  i1 = point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
  i2 = point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
  p = PolyArea([o, i1, i2])
  m = !p
  r = rings(m)
  @test m isa MultiPolygon
  @test length(r) == 4
  g = parent(m)
  @test length(g) == 3
  @test rings(g[1])[1] ≈ Ring(point.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test rings(g[1])[2] == Ring(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))
  @test rings(g[2]) == [Ring(point.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))]
  @test rings(g[3]) == [Ring(point.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))]

  b = Box(point(0, 0), point(1, 1))
  p = !b
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(point.([(0 - τ, 0 - τ), (1 + τ, 0 - τ), (1 + τ, 1 + τ), (0 - τ, 1 + τ)]))
  @test r[2] == Ring(point.([(0, 0), (0, 1), (1, 1), (1, 0)]))

  b = Ball(point(0, 0), T(1))
  p = !b
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] ≈ Ring(point.([(-1 - τ, -1 - τ), (1 + τ, -1 - τ), (1 + τ, 1 + τ), (-1 - τ, 1 + τ)]))
end
