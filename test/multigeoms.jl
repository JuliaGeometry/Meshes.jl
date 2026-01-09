@testitem "MultiGeometry" setup = [Setup] begin
  outer = cart.([(0, 0), (1, 0), (1, 1), (0, 1)])
  hole1 = cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)])
  hole2 = cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)])
  poly = PolyArea([outer, hole1, hole2])
  multi = Multi([poly, poly])
  @test multi == multi
  @test multi ≈ multi
  @test paramdim(multi) == 2
  @test crs(multi) <: Cartesian{NoDatum}
  @test Meshes.lentype(multi) == ℳ
  @test vertex(multi, 1) == vertex(poly, 1)
  @test vertices(multi) == [vertices(poly); vertices(poly)]
  @test nvertices(multi) == nvertices(poly) + nvertices(poly)
  @test rings(multi) == [rings(poly); rings(poly)]

  poly1 = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  poly2 = PolyArea(cart.([(1, 1), (2, 1), (2, 2), (1, 2)]))
  multi = Multi([poly1, poly2])
  @test vertices(multi) == [vertices(poly1); vertices(poly2)]
  @test nvertices(multi) == nvertices(poly1) + nvertices(poly2)
  @test area(multi) == area(poly1) + area(poly2)
  @test perimeter(multi) == perimeter(poly1) + perimeter(poly2)
  @test centroid(multi) == cart(1, 1)
  @test cart(0.5, 0.5) ∈ multi
  @test cart(1.5, 1.5) ∈ multi
  @test cart(1.5, 0.5) ∉ multi
  @test cart(0.5, 1.5) ∉ multi
  @test sprint(show, multi) == "Multi(2×PolyArea)"
  @test sprint(show, MIME"text/plain"(), multi) == """
  Multi(2×PolyArea)
  ├─ PolyArea((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))
  └─ PolyArea((x: 1.0 m, y: 1.0 m), ..., (x: 1.0 m, y: 2.0 m))"""

  box1 = Box(cart(0, 0), cart(1, 1))
  box2 = Box(cart(1, 1), cart(2, 2))
  mbox = Multi([box1, box2])
  @test sprint(show, mbox) == "Multi(2×Box)"
  @test sprint(show, MIME"text/plain"(), mbox) == """
  Multi(2×Box)
  ├─ Box(min: (x: 0.0 m, y: 0.0 m), max: (x: 1.0 m, y: 1.0 m))
  └─ Box(min: (x: 1.0 m, y: 1.0 m), max: (x: 2.0 m, y: 2.0 m))"""

  box1 = Box(cart(0, 0), cart(1, 1))
  box2 = Box(cart(1, 1), cart(2, 2))
  mbox = Multi([box1, box2])
  equaltest(mbox)
  isapproxtest(mbox)

  # constructor with iterator
  grid = cartgrid(10, 10)
  multi = Multi(grid)
  @test parent(multi) == collect(grid)

  # unique vertices
  poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  quad = Quadrangle(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  multi = Multi([poly, quad])
  @test unique(multi) == multi
  @test sprint(show, multi) == "Multi(2×Polygon)"
  @test sprint(show, MIME"text/plain"(), multi) == """
  Multi(2×Polygon)
  ├─ PolyArea((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))
  └─ Quadrangle((x: 0.0 m, y: 0.0 m), ..., (x: 0.0 m, y: 1.0 m))"""

  # type aliases
  p = cart(0, 0)
  segm = Segment(cart(0, 0), cart(1, 1))
  rope = Rope(cart.([(0, 0), (1, 0), (1, 1)]))
  ring = Ring(cart.([(0, 0), (1, 0), (1, 1)]))
  tri = Triangle(cart(0, 0), cart(1, 0), cart(1, 1))
  poly = PolyArea(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  @test Multi([p, p]) isa MultiPoint
  @test Multi([segm, segm]) isa MultiSegment
  @test Multi([rope, ring]) isa MultiChain
  @test Multi([rope, rope]) isa MultiRope
  @test Multi([ring, ring]) isa MultiRing
  @test Multi([tri, tri]) isa MultiPolygon
  @test Multi([poly, poly]) isa MultiPolygon

  # CRS propagation
  poly1 = PolyArea(merc.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  poly2 = PolyArea(merc.([(1, 1), (2, 1), (2, 2), (1, 2)]))
  multi = Multi([poly1, poly2])
  @test crs(centroid(multi)) === crs(multi)

  # vertex iteration
  ring1 = Ring(cart.([(0, 0), (1, 0), (1, 1), (0, 1)]))
  ring2 = Ring(cart.([(0, 0), (2, 0), (2, 2), (0, 2)]))
  ring3 = Ring(cart.([(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]))
  ring4 = Ring(cart.([(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]))
  poly1 = PolyArea(ring1)
  poly2 = PolyArea(ring2)
  poly3 = PolyArea([ring1, ring3])
  poly4 = PolyArea([ring2, ring4])
  multi1 = Multi([ring1, ring2, ring3, ring4])
  multi2 = Multi([poly1, poly2])
  multi3 = Multi([poly3, poly4])
  vertextest(multi1)
  vertextest(multi2)
  vertextest(multi3)
end
