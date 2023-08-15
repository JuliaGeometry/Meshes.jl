@testset "Clipping" begin
  window = Ring(Point.([(0,0), (0,4), (5,4),(5,0)]))
  poly = PolyArea(Point.([(0,2), (3,5), (6,2)]))

  clippedpoly = clip(poly, window, SutherlandHodgman())

  @test clippedpoly ≈ PolyArea(Point.([(0,2), (2,4), (4,4), (5,3), (5,2)]))
end
