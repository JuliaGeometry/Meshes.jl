@testset "PointSets" begin
  pset = PointSet(rand(P2, 100))
  @test embeddim(pset) == 2
  @test coordtype(pset) == T

  pset = PointSet(rand(P3, 100))
  @test embeddim(pset) == 3
  @test coordtype(pset) == T

  pset1 = PointSet(P3(1,2,3), P3(4,5,6))
  pset2 = PointSet([T.((1,2,3)), T.((4,5,6))])
  pset3 = PointSet(T.((1,2,3)), T.((4,5,6)))
  @test pset1 == pset2 == pset3
end
