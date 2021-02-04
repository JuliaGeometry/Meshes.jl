@testset "PointSets" begin
  pset = PointSet(rand(P2, 100))
  @test embeddim(pset) == 2
  @test coordtype(pset) == T

  pset = PointSet(rand(P3, 100))
  @test embeddim(pset) == 3
  @test coordtype(pset) == T
end
