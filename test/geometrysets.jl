@testset "GeometrySet" begin
  @test !isgrid(GeometrySet)

  s = Segment(P2(0,0), P2(1,1))
  t = Triangle(P2(0,0), P2(1,0), P2(0,1))
  p = PolyArea(P2[(0,0), (1,0), (1,1), (0,1), (0,0)])
  gset = GeometrySet([s, t, p])
  @test [centroid(gset, i) for i in 1:3] == P2[(1/2,1/2), (1/3,1/3), (1/2,1/2)]
end
