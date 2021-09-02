@testset "Convex Hull" begin
  
  # Upper part of the Convex Hull
  pset1 = PointSet(P2(0,0),P2(1,1),P2(1.5,1.0))
  pset2 = PointSet(P2(0,0),P2(1,1),P2(1.5,1.0),P2(1.5,3.))
  pset3 = PointSet(P2(0,0),P2(1,1),P2(1.5,1.0),P2(1.5,3.),P2(2,3))
  pset4 = PointSet(P2(0,0),P2(1,1),P2(1.5,1.0),P2(1.5,3.),P2(2,3),P2(2,4))

  @test upperhull(pset1) == [[0.0, 0.0],[1.0, 1.0],[1.5, 1.0]]
  @test upperhull(pset2) == [[0.0, 0.0],[1.5, 3.0]]
  @test upperhull(pset3) == [[0.0, 0.0],[1.5, 3.0],[2.0, 3.0]]
  @test upperhull(pset4) == [[0.0, 0.0],[1.5, 3.0],[2.0, 4.0]]

end