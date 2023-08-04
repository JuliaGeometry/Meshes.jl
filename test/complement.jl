@testset "complement" begin
  t = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
  p = !t
  @test p isa PolyArea
  @test rings(p) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)]), Ring(P2[(0, 0), (1, 1), (1, 0)])]

  q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
  p = !q
  @test p isa PolyArea
  @test rings(p) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)]), Ring(P2[(0, 0), (0, 1), (1, 1), (1, 0)])]

  p = PolyArea(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
  n = !p
  @test n isa PolyArea
  @test rings(n) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)]), Ring(P2[(0, 0), (0, 1), (1, 1), (1, 0)])]

  o = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
  i1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
  i2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
  p = PolyArea(o, [i1, i2])
  m = !p
  @test m isa MultiPolygon
  g = collect(m)
  @test length(g) == 3
  @test rings(g[1]) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)]), Ring(P2[(0, 0), (0, 1), (1, 1), (1, 0)])]
  @test rings(g[2]) == [Ring(P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])]
  @test rings(g[3]) == [Ring(P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])]

  b = Box(P2(0, 0), P2(1, 1))
  p = !b
  @test p isa PolyArea
  rings(p)
  @test rings(p) == [Ring(P2[(0, 0), (1, 0), (1, 1), (0, 1)]), Ring(P2[(0, 0), (0, 1), (1, 1), (1, 0)])]

  b = Ball(P2(0, 0), T(1))
  p = !b
  r = rings(p)
  @test p isa PolyArea
  @test length(r) == 2
  @test r[1] == Ring(P2[(-1, -1), (1, -1), (1, 1), (-1, 1)])
end
