@testset "Multi" begin
  outer = P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
  hole1 = P2[(0.2, 0.2), (0.4, 0.2), (0.4, 0.4), (0.2, 0.4), (0.2, 0.2)]
  hole2 = P2[(0.6, 0.2), (0.8, 0.2), (0.8, 0.4), (0.6, 0.4), (0.6, 0.2)]
  poly = PolyArea(outer, [hole1, hole2])
  multi = Multi([poly, poly])
  @test multi == multi
  @test paramdim(multi) == 2
  @test vertex(multi, 1) == vertex(poly, 1)
  @test vertices(multi) == [vertices(poly); vertices(poly)]
  @test nvertices(multi) == nvertices(poly) + nvertices(poly)
  @test boundary(multi) == merge(boundary(poly), boundary(poly))
  @test chains(multi) == [chains(poly); chains(poly)]

  poly1 = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
  poly2 = PolyArea(P2[(1, 1), (2, 1), (2, 2), (1, 2), (1, 1)])
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

  box1 = Box(P2(0, 0), P2(1, 1))
  box2 = Box(P2(1, 1), P2(2, 2))
  mbox = Multi([box1, box2])
  mchn = boundary(mbox)
  noth = boundary(mchn)
  @test mchn isa Multi
  @test isnothing(noth)
  @test length(mchn) == T(8)

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
end
