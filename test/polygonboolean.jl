@testitem "Polygon Boolean internals" setup = [Setup] begin
  # _insertintersections! function test
  # simple case: two rings intersecting at two points
  ring1 = Ring(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  ring2 = Ring(cart(2, -1), cart(2, 5), cart(5, 5), cart(5, -1))
  rings = [ring1, ring2]

  segs = segments.(rings)
  allsegs = Iterators.flatten(segs)
  intersections, seginds = Meshes.pairwiseintersect(allsegs)

  Meshes._insertintersections!(intersections, seginds, vertices.(rings))

  # check that intersections were inserted correctly
  newring1 = rings[1]
  newring2 = rings[2]

  @test length(vertices(newring1)) == 6  # two intersections added
  @test length(vertices(newring2)) == 6  # two intersections added

  @test cart(2, 0) in vertices(newring1)
  @test cart(2, 4) in vertices(newring1)
  @test cart(2, 0) in vertices(newring2)
  @test cart(2, 4) in vertices(newring2)

  # degenerate case: intersections at segment endpoints
  ring3 = Ring(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  ring4 = Ring(cart(4, 0), cart(4, 4), cart(8, 4), cart(8, 0))
  rings2 = [ring3, ring4]
  segs2 = segments.(rings2)
  allsegs2 = Iterators.flatten(segs2)
  intersections2, seginds2 = Meshes.pairwiseintersect(allsegs2)
  Meshes._insertintersections!(intersections2, seginds2, vertices.(rings2))
  newring3 = rings2[1]
  newring4 = rings2[2]
  @test length(vertices(newring3)) == 4  # no new points added
  @test length(vertices(newring4)) == 4  # no new points added

  ## _filled function test
  # test fills above
  @test Meshes._filled(union, 0b0001, true) == true
  @test Meshes._filled(union, 0b0010, true) == false
  @test Meshes._filled(union, 0b0000, true) == false
  @test Meshes._filled(intersect, 0b0001, true) == false
  @test Meshes._filled(intersect, 0b0100, true) == false
  @test Meshes._filled(intersect, 0b0101, true) == true
  @test Meshes._filled(setdiff, 0b0001, true) == true
  @test Meshes._filled(setdiff, 0b0101, true) == false
  @test Meshes._filled(symdiff, 0b0001, true) == true
  @test Meshes._filled(symdiff, 0b0100, true) == true
  @test Meshes._filled(symdiff, 0b0110, true) == true
  # test fills below
  @test Meshes._filled(union, 0b0100, false) == false
  @test Meshes._filled(union, 0b1000, false) == true
  @test Meshes._filled(union, 0b0000, false) == false
  @test Meshes._filled(intersect, 0b0010, false) == false
  @test Meshes._filled(intersect, 0b1010, false) == true
  @test Meshes._filled(setdiff, 0b0010, false) == true
  @test Meshes._filled(setdiff, 0b1010, false) == false
  @test Meshes._filled(symdiff, 0b1000, false) == true
  @test Meshes._filled(symdiff, 0b0110, false) == true
  @test Meshes._filled(symdiff, 0b1110, false) == false
end
