@testset "Bounding boxes" begin
  @test boundingbox((0,0),1) == ((-1,-1),(1,1))
  @test boundingbox((1,1),1) == ((0,0),(2,2))
end
