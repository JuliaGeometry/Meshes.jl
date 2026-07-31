@testitem "segments" setup = [Setup] begin
  # segments of chains
  c = Rope(cart.([(1, 1), (2, 2), (3, 3)]))
  @test collect(segments(c)) == [Segment(cart(1, 1), cart(2, 2)), Segment(cart(2, 2), cart(3, 3))]
  c = Ring(cart.([(1, 1), (2, 2), (3, 3)]))
  @test collect(segments(c)) ==
        [Segment(cart(1, 1), cart(2, 2)), Segment(cart(2, 2), cart(3, 3)), Segment(cart(3, 3), cart(1, 1))]

  # degenerate rings
  r = Ring(cart.([(0, 0)]))
  @test collect(segments(r)) == [Segment(cart(0, 0), cart(0, 0))]
  r = Ring(cart.([(0, 0), (1, 1)]))
  @test collect(segments(r)) == [Segment(cart(0, 0), cart(1, 1)), Segment(cart(1, 1), cart(0, 0))]

  # segments of meshes
  g = cartgrid(2)
  @test length(collect(segments(g))) == 2
  g = cartgrid(2, 2)
  h = topoconvert(HalfEdgeTopology, g)
  s = topoconvert(SimpleTopology, g)
  @test length(collect(segments(g))) == 12
  @test length(collect(segments(h))) == 12
  @test length(collect(segments(s))) == 12
  g = cartgrid(2, 3)
  h = topoconvert(HalfEdgeTopology, g)
  s = topoconvert(SimpleTopology, g)
  @test length(collect(segments(g))) == 17
  @test length(collect(segments(h))) == 17
  @test length(collect(segments(s))) == 17
end
