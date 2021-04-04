@testset "Distances" begin
  p = P2(0, 1)
  l = Line(P2(0, 0), P2(1, 0))
  @test evaluate(Euclidean(), p, l) == T(1)
  @test evaluate(Euclidean(), l, p) == T(1)

  p1, p2 = P2(1, 0), P2(0, 1)
  @test Chebyshev()(p1, p2) == T(1)
  @test evaluate(Chebyshev(), p1, p2) == T(1)

  p = P2(68, 259)
  l = Line(P2(68, 260), P2(69, 261))
  @test evaluate(Euclidean(), p, l) ≤ T(0.8)
end
