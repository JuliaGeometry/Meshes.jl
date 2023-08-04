@testset "Multi" begin
  outer = P2[(0, 0), (1, 0), (1, 1), (0, 1)]
  hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4)]
  hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4)]
  poly = PolyArea(outer, [hole1, hole2])
  multi = Multi([poly, poly])
  @test multi == multi
  @test multi ≈ multi
  @test paramdim(multi) == 2
  @test vertex(multi, 1) == vertex(poly, 1)
  @test vertices(multi) == [vertices(poly); vertices(poly)]
  @test nvertices(multi) == nvertices(poly) + nvertices(poly)
  @test boundary(multi) == merge(boundary(poly), boundary(poly))
  @test rings(multi) == [rings(poly); rings(poly)]

  poly1 = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  poly2 = PolyArea(P2[(1, 1), (2, 1), (2, 2), (1, 2)])
  multi = Multi([poly1, poly2])
  @test vertices(multi) == [vertices(poly1); vertices(poly2)]
  @test nvertices(multi) == nvertices(poly1) + nvertices(poly2)
  @test area(multi) == area(poly1) + area(poly2)
  @test perimeter(multi) == perimeter(poly1) + perimeter(poly2)
  @test centroid(multi) == P2(1, 1)
  @test P2(0.5, 0.5) ∈ multi
  @test P2(1.5, 1.5) ∈ multi
  @test P2(1.5, 0.5) ∉ multi
  @test P2(0.5, 1.5) ∉ multi
  if T == Float32
    @test sprint(show, multi) == "MultiPolyArea{2,Float32}"
    @test sprint(show, MIME"text/plain"(), multi) ==
          "MultiPolyArea{2,Float32}\n  └─PolyArea(4-Ring)\n  └─PolyArea(4-Ring)"
  elseif T == Float64
    @test sprint(show, multi) == "MultiPolyArea{2,Float64}"
    @test sprint(show, MIME"text/plain"(), multi) ==
          "MultiPolyArea{2,Float64}\n  └─PolyArea(4-Ring)\n  └─PolyArea(4-Ring)"
  end

  box1 = Box(P2(0, 0), P2(1, 1))
  box2 = Box(P2(1, 1), P2(2, 2))
  mbox = Multi([box1, box2])
  mchn = boundary(mbox)
  noth = boundary(mchn)
  @test mchn isa Multi
  @test isnothing(noth)
  @test length(mchn) == T(8)
  if T == Float32
    @test sprint(show, mbox) == "MultiBox{2,Float32}"
    @test sprint(show, MIME"text/plain"(), mbox) ==
          "MultiBox{2,Float32}\n  └─Box{2, Float32}(Point(0.0f0, 0.0f0), Point(1.0f0, 1.0f0))\n  └─Box{2, Float32}(Point(1.0f0, 1.0f0), Point(2.0f0, 2.0f0))"
  elseif T == Float64
    @test sprint(show, mbox) == "MultiBox{2,Float64}"
    @test sprint(show, MIME"text/plain"(), mbox) ==
          "MultiBox{2,Float64}\n  └─Box{2, Float64}(Point(0.0, 0.0), Point(1.0, 1.0))\n  └─Box{2, Float64}(Point(1.0, 1.0), Point(2.0, 2.0))"
  end

  # constructor with iterator
  grid = CartesianGrid{T}(10, 10)
  multi = Multi(grid)
  @test collect(multi) == collect(grid)

  # boundary of multi-3D-geometry
  box1 = Box(P3(0, 0, 0), P3(1, 1, 1))
  box2 = Box(P3(1, 1, 1), P3(2, 2, 2))
  mbox = Multi([box1, box2])
  mesh = boundary(mbox)
  @test mesh isa Mesh
  @test nvertices(mesh) == 16
  @test nelements(mesh) == 12

  # unique vertices
  poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  quad = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
  multi = Multi([poly, quad])
  @test unique(multi) == multi
  if T == Float32
    @test sprint(show, multi) == "MultiPolygon{2,Float32}"
    @test sprint(show, MIME"text/plain"(), multi) ==
          "MultiPolygon{2,Float32}\n  └─PolyArea(4-Ring)\n  └─Quadrangle(Point(0.0f0, 0.0f0), Point(1.0f0, 0.0f0), Point(1.0f0, 1.0f0), Point(0.0f0, 1.0f0))"
  elseif T == Float64
    @test sprint(show, multi) == "MultiPolygon{2,Float64}"
    @test sprint(show, MIME"text/plain"(), multi) ==
          "MultiPolygon{2,Float64}\n  └─PolyArea(4-Ring)\n  └─Quadrangle(Point(0.0, 0.0), Point(1.0, 0.0), Point(1.0, 1.0), Point(0.0, 1.0))"
  end

  # type aliases
  point = P2(0, 0)
  segm = Segment(P2(0, 0), P2(1, 1))
  rope = Rope(P2[(0, 0), (1, 0), (1, 1)])
  ring = Ring(P2[(0, 0), (1, 0), (1, 1)])
  tri = Triangle(P2(0, 0), P2(1, 0), P2(1, 1))
  poly = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
  @test Multi([point, point]) isa MultiPoint
  @test Multi([segm, segm]) isa MultiSegment
  @test Multi([rope, rope]) isa MultiRope
  @test Multi([ring, ring]) isa MultiRing
  @test Multi([tri, tri]) isa MultiPolygon
  @test Multi([poly, poly]) isa MultiPolygon
end
