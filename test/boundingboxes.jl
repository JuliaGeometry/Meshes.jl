@testset "Bounding boxes" begin
  # sphere
  @test boundingbox(Sphere(P2(0,0),T(1))) == Box(P2(-1,-1),P2(1,1))
  @test boundingbox(Sphere(P2(1,1),T(1))) == Box(P2(0,0),P2(2,2))

  # point set
  @test boundingbox(PointSet(T[0 1 2; 0 2 1])) == Box(P2(0,0), P2(2,2))
  @test boundingbox(PointSet(T[1 2; 2 1])) == Box(P2(1,1), P2(2,2))

  # cartesian grid
  @test boundingbox(CartesianGrid{T}(10,10)) == Box(P2(0,0), P2(10,10))
  @test boundingbox(CartesianGrid{T}(100,200)) == Box(P2(0,0), P2(100,200))
  @test boundingbox(CartesianGrid((10,10), T.((1,1)), T.((1,1)))) == Box(P2(1,1), P2(11,11))
end
