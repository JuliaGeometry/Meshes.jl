@testset "Multi" begin
  outer = P2[(0,0),(1,0),(1,1),(0,1),(0,0)]
  hole1 = P2[(0.2,0.2),(0.4,0.2),(0.4,0.4),(0.2,0.4),(0.2,0.2)]
  hole2 = P2[(0.6,0.2),(0.8,0.2),(0.8,0.4),(0.6,0.4),(0.6,0.2)]
  poly  = PolyArea(outer, [hole1, hole2])
  multi = Multi([poly, poly])
  @test multi == multi
  @test chains(multi) == [chains(poly); chains(poly)]
end
