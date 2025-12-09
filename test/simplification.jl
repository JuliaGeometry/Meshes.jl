@testitem "Selinger" setup = [Setup] begin
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)]))
  s1 = simplify(c, SelingerSimplification(T(0.1)))
  s2 = simplify(c, SelingerSimplification(T(0.5)))
  @test s1 == Ring(cart.([(1, 0), (1, 1), (2, 1), (2, 2), (0, 2), (0, 0)]))
  @test s2 == Ring(cart.([(1, 0), (2, 2), (0, 2), (0, 0)]))

  b = Box(cart(0, 0), cart(1, 1))
  s = simplify(b, SelingerSimplification(T(0.1)))
  @test nvertices(s) == 4
  s = simplify(b, SelingerSimplification(T(0.8)))
  @test nvertices(s) == 3
end

@testitem "DouglasPeucker" setup = [Setup] begin
  c = Ring(cart.([(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)]))
  s1 = simplify(c, DouglasPeuckerSimplification(T(0.1)))
  s2 = simplify(c, DouglasPeuckerSimplification(T(0.5)))
  @test s1 == Ring(cart.([(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)]))
  @test s2 == Ring(cart.([(0, 0), (1.5, 0.5), (0, 1)]))

  p = PolyArea(Ring(cart.([(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)])))
  s1 = simplify(p, DouglasPeuckerSimplification(T(0.5)))
  @test s1 == PolyArea(Ring(cart.([(0, 0), (1.5, 0.5), (0, 1)])))
  m = Multi([p, p])
  s2 = simplify(m, DouglasPeuckerSimplification(T(0.5)))
  @test s2 == Multi([s1, s1])
  d = GeometrySet([p, p])
  s3 = simplify(d, DouglasPeuckerSimplification(T(0.5)))
  @test s3 == GeometrySet([s1, s1])
end

@testitem "MinMax" setup = [Setup] begin
  # Selinger
  c = Ring(cart.([(0, 0), (1, 0), (1, 1), (2, 1), (2, 2), (1, 2), (0, 2), (0, 1)]))
  s1 = simplify(c, SelingerSimplification(T(0.1)))
  s2 = simplify(c, MinMaxSimplification(SelingerSimplification, max=6))
  @test nvertices(s2) ≤ nvertices(s1)
  s1 = simplify(c, SelingerSimplification(T(0.5)))
  s2 = simplify(c, MinMaxSimplification(SelingerSimplification, max=4))
  @test nvertices(s2) ≤ nvertices(s1)

  # Douglas-Peucker
  c = Ring(cart.([(0, 0), (1, 0), (1.5, 0.5), (1, 1), (0, 1)]))
  s1 = simplify(c, DouglasPeuckerSimplification(T(0.1)))
  s2 = simplify(c, MinMaxSimplification(DouglasPeuckerSimplification, max=6))
  @test s1 ≗ s2
  s1 = simplify(c, DouglasPeuckerSimplification(T(0.5)))
  s2 = simplify(c, MinMaxSimplification(DouglasPeuckerSimplification, max=4))
  @test s1 ≗ s2
end
