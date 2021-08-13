@testset "upperhull" begin
  
  pset1=PointSet((0.,0.),(1.,1.),(1.5,1.0))
  
  @test upperhull(pset1) = [[0.0, 0.0],[1.0, 1.0], [1.5, 1.0]]
  
end