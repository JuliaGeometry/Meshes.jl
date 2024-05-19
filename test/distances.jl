@testset "Distances" begin
  p = point(0, 1)
  l = Line(point(0, 0), point(1, 0))
  @test evaluate(Euclidean(), p, l) == T(1) * u"m"
  @test evaluate(Euclidean(), l, p) == T(1) * u"m"

  p1, p2 = point(1, 0), point(0, 1)
  @test evaluate(Chebyshev(), p1, p2) == T(1) * u"m"

  p = point(68, 259)
  l = Line(point(68, 260), point(69, 261))
  @test evaluate(Euclidean(), p, l) ≤ T(0.8) * u"m"

  line1 = Line(point(-1, 0, 0), point(1, 0, 0))
  line2 = Line(point(0, -1, 1), point(0, 1, 1))  # line2 ⟂ line1, z++
  line3 = Line(point(-1, 1, 0), point(1, 1, 0))  # line3 ∥ line1
  line4 = Line(point(-2, 0, 0), point(2, 0, 0))  # line4 colinear with line1
  line5 = Line(point(0, -1, 0), point(0, 1, 0))  # line5 intersects line1
  @test evaluate(Euclidean(), line1, line2) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line3) ≈ T(1) * u"m"
  @test evaluate(Euclidean(), line1, line4) ≈ T(0) * u"m"
  @test evaluate(Euclidean(), line1, line5) ≈ T(0) * u"m"
end
