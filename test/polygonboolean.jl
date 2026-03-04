@testitem "Polygon Boolean internals" setup = [Setup] begin
  # insertintersections! function test
  # simple case: two rings intersecting at two points
  ring1 = Ring(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  ring2 = Ring(cart(2, -1), cart(2, 5), cart(5, 5), cart(5, -1))
  rings = [ring1, ring2]

  segs = segments.(rings)
  allsegs = Iterators.flatten(segs)
  intersections, seginds = Meshes.pairwiseintersect(allsegs)

  Meshes._insertintersections!(rings, intersections, seginds)

  @test length(vertices(rings[1])) == 6  # two intersections added
  @test length(vertices(rings[2])) == 6  # two intersections added

  @test all(vertices(rings[1]) .== cart.([(0, 0), (2, 0), (4, 0), (4, 4), (2, 4), (0, 4)]))
  @test all(vertices(rings[2]) .== cart.([(2, -1), (2, 0), (2, 4), (2, 5), (5, 5), (5, -1)]))

  # degenerate case: intersections at segment endpoints
  ring3 = Ring(cart(0, 0), cart(4, 0), cart(4, 4), cart(0, 4))
  ring4 = Ring(cart(4, 0), cart(4, 4), cart(8, 4), cart(8, 0))
  rings2 = [ring3, ring4]
  segs2 = segments.(rings2)
  allsegs2 = Iterators.flatten(segs2)
  intersections2, seginds2 = Meshes.pairwiseintersect(allsegs2)
  Meshes._insertintersections!(rings2, intersections2, seginds2)
  @test length(vertices(rings2[1])) == 4  # no new points added
  @test length(vertices(rings2[2])) == 4  # no new points added

  # _filledabove and _filledbelow correctness across all operations and bit patterns
  ops = [union, intersect, setdiff, symdiff]

  for op in ops
    for bits in 0b0000:0b1111
      result_above = Meshes._filledabove(op, bits)
      result_below = Meshes._filledbelow(op, bits)

      @test isa(result_above, Bool)
      @test isa(result_below, Bool)

      if op == union
        babove = (bits & 0b0001) == 1 || (bits & 0b0100) == 4
        bbelow = (bits & 0b0010) == 2 || (bits & 0b1000) == 8
      elseif op == intersect
        babove = (bits & 0b0001) == 1 && (bits & 0b0100) == 4
        bbelow = (bits & 0b0010) == 2 && (bits & 0b1000) == 8
      elseif op == setdiff
        babove = (bits & 0b0001) == 1 && (bits & 0b0100) != 4
        bbelow = (bits & 0b0010) == 2 && (bits & 0b1000) != 8
      elseif op == symdiff
        babove = ((bits & 0b0001) == 1) ⊻ ((bits & 0b0100) == 4)
        bbelow = ((bits & 0b0010) == 2) ⊻ ((bits & 0b1000) == 8)
      end

      @test result_above == babove
      @test result_below == bbelow
    end
  end
end
