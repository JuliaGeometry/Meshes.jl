@testset "Bounding boxes" begin
  @test boundingbox(Sphere(P2(0,0),T(1))) == Box(P2(-1,-1),P2(1,1))
  @test boundingbox(Sphere(P2(1,1),T(1))) == Box(P2(0,0),P2(2,2))
  @test boundingbox(CartesianGrid{T}(10,10)) == Box(P2(0,0), P2(10,10))
end
