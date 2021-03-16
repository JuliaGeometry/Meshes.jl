@testset "Distances" begin
  p = P2(0,1)
  l = Line(P2(0,0), P2(1,0))
  @test evaluate(Euclidean(), p, l) == T(1)
  @test evaluate(Euclidean(), l, p) == T(1)
end
