@testset "Bounding boxes" begin
  @test boundingbox(Sphere((0,0),1)) == Box(Point(-1,-1),Point(1,1))
  @test boundingbox(Sphere((1,1),1)) == Box(Point(0,0),Point(2,2))
end
