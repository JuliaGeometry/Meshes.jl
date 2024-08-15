@testset "TransformedGeometry" begin
  b = Box(cart(0, 0), cart(1, 1))
  t = Translate(T(1), T(2))
  tb = Meshes.TransformedGeometry(b, t)
  @test parent(tb) == b
  @test Meshes.transform(tb) == t
  t2 = Scale(T(2), T(3))
  tb2 = Meshes.TransformedGeometry(tb, t2)
  @test Meshes.transform(tb2) == (t → t2)
  @test paramdim(tb) == paramdim(b)
  @test centroid(tb) == t(centroid(b))
  @test discretize(tb) == t(discretize(b))
  @test tb(T(0.5), T(0.5)) == t(b(T(0.5), T(0.5)))
  t3 = Scale(T(2), T(2))
  tb3 = Meshes.TransformedGeometry(b, t3)
  @test measure(tb3) == 4 * measure(b)
  equaltest(tb)
  isapproxtest(tb)

  b = Ball(latlon(0, 0), T(1))
  t = Proj(Cartesian)
  tb = Meshes.TransformedGeometry(b, t)
  @test paramdim(tb) == paramdim(b)
  @test centroid(tb) == t(centroid(b))

  s = Sphere(latlon(0, 0), T(1))
  t = Proj(Cartesian)
  ts = Meshes.TransformedGeometry(s, t)
  @test paramdim(ts) == paramdim(s)
  @test centroid(ts) == t(centroid(s))

  s = Segment(cart(0, 0), cart(1, 1))
  t = Translate(T(1), T(2))
  ts = Meshes.TransformedGeometry(s, t)
  @test vertex(ts, 1) == t(vertex(s, 1))
  @test vertices(ts) == t.(vertices(s))
  @test nvertices(ts) == nvertices(s)

  p = PolyArea(cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  t = Translate(T(1), T(2))
  tp = Meshes.TransformedGeometry(p, t)
  @test vertex(tp, 1) == t(vertex(p, 1))
  @test vertices(tp) == t.(vertices(p))
  @test nvertices(tp) == nvertices(p)
  @test rings(tp) == t.(rings(p))
  p2 = PolyArea(cart(0, 0), cart(0, 0), cart(1, 0), cart(1, 1), cart(0, 1))
  tp2 = Meshes.TransformedGeometry(p2, t)
  @test unique(tp2) == tp

  b = Box(cart(0, 0), cart(1, 1))
  t = Translate(T(1), T(2))
  tb = Meshes.TransformedGeometry(b, t)
  @test sprint(show, tb) ==
        "TransformedBox(geometry: Box(min: (x: 0.0 m, y: 0.0 m), max: (x: 1.0 m, y: 1.0 m)), transform: Translate(offsets: (1.0 m, 2.0 m)))"
  @test sprint(show, MIME"text/plain"(), tb) == """
  TransformedBox
  ├─ geometry: Box(min: (x: 0.0 m, y: 0.0 m), max: (x: 1.0 m, y: 1.0 m))
  └─ transform: Translate(offsets: (1.0 m, 2.0 m))"""
end
